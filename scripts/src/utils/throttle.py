import logging
import random
import sqlite3
import time
from pprint import pformat

from .config import Config
from .errors import APIError

config = Config()
logger = logging.getLogger(__name__)


class ENDPOINTS:
    GBIF_SLOW = {
        'requests_per_second': 1,
        'name': 'gbif_slow',
    }
    GBIF_FAST = {
        'requests_per_second': 10,
        'name': 'gbif_fast',
    }
    ENTREZ = {
        'requests_per_second': 10,
        'name': 'entrez',
    }
    BOLD = {
        'requests_per_second': 10,
        'name': 'bold',
    }


class Throttle:
    """Use SQLite2 database to coordinate throttling of API requests.

    This is necessary to avoid hitting API rate limits, or overwhelming the
    server. Each endpoint (identified by name) is throttled independently to
    allow for request rates to be set per-service, and to for throttles to be
    managed independently.

    The endpoint arg should be a dict of:
        {
          'requests_per_second': int,  # Max requests per second
          'name': str,                 # Name to identify this endpoint
        }
    """

    FIELD_NAME = 'timestamp'

    def __init__(
        self,
        endpoint: dict,
    ):
        self.requests_per_second = endpoint['requests_per_second']
        self.db_path = config.throttle_sqlite_path
        self.name = endpoint['name']
        self._initialize_db()

    def __enter__(self):
        self.await_release()

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def _initialize_db(self):
        """Create table for tracking request timestamps."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(f"""
                CREATE TABLE IF NOT EXISTS {self.name} (
                    {self.FIELD_NAME} INTEGER
                )
            """)
            conn.commit()

    def await_release(self):
        """Query sqlite DB for permission to send a request.

        The DB table keeps track of requests sent across processes by writing
        a timestamp for each request. If the number of requests in the last
        second exceeds the limit, the request is blocked until the sliding
        window is clear.
        """
        while True:
            with sqlite3.connect(self.db_path, isolation_level=None) as conn:
                try:
                    # Lock the database for writing
                    conn.execute("BEGIN IMMEDIATE")

                    now = int(time.time() * 1000)
                    window_start = now - 2000  # Two-second sliding window

                    # Remove expired timestamps (older than 1s)
                    conn.execute(
                        f"DELETE FROM {self.name} WHERE {self.FIELD_NAME} < ?",
                        (window_start,))

                    # Count remaining requests in the last second
                    request_count = conn.execute(
                        f"SELECT COUNT(*) FROM {self.name}"
                    ).fetchone()[0]

                    if request_count < self.requests_per_second:
                        # Insert current timestamp atomically
                        conn.execute(
                            f"INSERT INTO {self.name} ({self.FIELD_NAME})"
                            " VALUES (?)",
                            (now,)
                        )
                        conn.commit()
                        return

                    # Rollback if the request limit is exceeded
                    conn.rollback()

                except sqlite3.OperationalError:
                    # Handle potential lock contention gracefully
                    pass

            # Sleep for a random interval to reduce race conditions collisions
            time.sleep(round(random.uniform(0.01, 0.1), 3))

    def with_retry(self, func, args=[], kwargs={}):
        retries = config.MAX_API_RETRIES
        while True:
            try:
                with self:
                    logger.debug("Throttle released. Sending request to"
                                 f" {func.__name__}...")
                    return func(*args, **kwargs)
            except Exception as exc:
                sleep_seconds = 1
                retries -= 1
                if '429' in str(exc):
                    sleep_seconds = 600
                    logger.warning(
                        "API rate limit exceeded. Waiting 10 minutes before"
                        " next retry.")
                    retries = config.MAX_API_RETRIES
                elif retries <= 0:
                    raise APIError(
                        'Failed to fetch data from API after'
                        f' {config.MAX_API_RETRIES} retries. Please try'
                        f' resuming this job at a later time.'
                        f'\nException: {exc}'
                    )
                logger.warning(
                    "Exception encountered in call to endpoint"
                    f" {func.__name__} Retrying {retries} more times."
                    f" Exception: {exc}\n"
                    f" Args:\n{pformat(args)}")
                time.sleep(sleep_seconds)
