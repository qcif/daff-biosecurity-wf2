import fcntl
import logging
import time
from pathlib import Path
from pprint import pformat

from .config import Config
from .errors import APIError

config = Config()
logger = logging.getLogger(__name__)


class Throttle:
    """Use filelock to throttle requests to the GBIF API.

    This is necessary to avoid hitting the rate limit of the GBIF API, or
    overwhelming the server.

    A custom interval and lock file can be set to allow for different
    throttling intervals for each endpoint.
    """

    def __init__(
        self,
        interval_sec: int,
        lock_file: Path,
    ):
        self.interval_sec = interval_sec
        self.lockfile = lock_file
        self._write_lockfile()

    def __enter__(self):
        self.acquire()

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def _write_lockfile(self):
        """Ensures the lock file exists."""
        if not self.lockfile.exists():
            with open(self.lockfile, "w") as f:
                f.write("")

    def acquire(self):
        while True:
            try:
                with open(self.lockfile, "r+") as f:
                    fcntl.flock(f, fcntl.LOCK_EX)  # Lock the file
                    time.sleep(self.interval_sec)
                    fcntl.flock(f, fcntl.LOCK_UN)  # Unlock the file
                    return
            except (FileNotFoundError, ValueError):
                self._write_lockfile()

            time.sleep(0.1)  # Wait before retrying

    def with_retry(self, func, args=[], kwargs={}):
        retries = config.MAX_API_RETRIES
        while True:
            try:
                with self:
                    logger.debug("Lock acquired. Sending request to"
                                 f" {func.__name__}...")
                    return func(*args, **kwargs)
            except Exception as exc:
                sleep_seconds = 1
                retries -= 1
                if '429' in str(exc):
                    sleep_seconds = 600
                    logger.warning(
                        "Entrez rate limit exceeded. Waiting 10 minutes before"
                        " next retry.")
                    retries = config.ENTREZ_MAX_RETRIES
                elif retries <= 0:
                    raise APIError(
                        'Failed to fetch data from GBIF API after'
                        f' {config.GBIF_MAX_RETRIES} retries. Please try'
                        f' resuming this job at a later time.'
                        f'\nException: {exc}'
                    )
                logger.warning(
                    "Exception encountered in call to GBIF endpoint"
                    f" {func.__name__} Retrying {retries} more times."
                    f" Exception: {exc}\n"
                    f" Args:\n{pformat(args)}")
                time.sleep(sleep_seconds)
