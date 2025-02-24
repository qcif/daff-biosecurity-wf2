"""Define error handling logic."""

from .config import Config

config = Config()

ERROR_CSV_HEADER = [
    "location",
    "message",
    "exception",
]


class APIError(Exception):
    """Raise this exception when an API request fails."""
    pass


class LOCATIONS:
    """These locations define where a reportable exception occurred and
    can be used to display it in the appropriate location in the report.
    """
    DATABASE_COVERAGE = 5.0
    DB_COVERAGE_TARGET = 5.1
    DB_COVERAGE_RELATED = 5.2
    DB_COVERAGE_RELATED_COUNTRY = 5.3


def write(location, msg, exc, query_dir=None):
    path = (query_dir or config.output_dir) / config.ERROR_FILENAME
    if not path.exists():
        path.write_text(','.join(ERROR_CSV_HEADER) + '\n')
    with path.open('a') as f:
        f.write(','.join[
            str(location),
            msg,
            str(exc),
        ] + '\n')
