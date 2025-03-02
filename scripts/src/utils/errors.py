"""Define error handling logic."""

import json

from .config import Config

config = Config()

ERROR_CSV_HEADER = [
    "location",
    "message",
    "exception",
    "data",
]


class APIError(Exception):
    """Raise this exception when an API request fails."""
    pass


class LOCATIONS:
    """These locations define where a reportable exception occurred and
    can be used to display it in the appropriate location in the report.
    """
    DATABASE_COVERAGE = 5.0
    DATABASE_COVERAGE_NO_GBIF_RECORD = 5.01
    DATABASE_COVERAGE_NO_TAXID = 5.02
    DB_COVERAGE_TARGET = 5.1
    DB_COVERAGE_RELATED = 5.2
    DB_COVERAGE_RELATED_COUNTRY = 5.3


def write(location, msg, exc, query_dir=None, data=None):
    parent = query_dir or config.output_dir
    next_path = parent / config.ERRORS_DIR / 'next.txt'
    next_path.parent.mkdir(parents=True, exist_ok=True)
    if next_path.exists():
        i = int(next_path.read_text())
    else:
        i = 1
    next_path.write_text(str(i + 1))
    path = parent / config.ERRORS_DIR / f'{i}.json'
    with path.open('w') as f:
        json.dump({
            "location": location,
            "message": msg,
            "exception": str(exc),
            "data": data,
        }, f)
