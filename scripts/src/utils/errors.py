"""Define error handling logic."""

import json

from .config import Config

config = Config()


class FASTAFormatError(Exception):
    """User provided a FASTA file that does not meet specifications."""

    def __init__(self, msg):
        super().__init__(f"FASTA format error: {msg}")


class MetadataFormatError(Exception):
    """User provided a metadata CSV file that does not meet specifications."""

    def __init__(self, msg):
        super().__init__(f"Metadata CSV format error: {msg}")


class APIError(Exception):
    """Raise this exception when an API request fails."""
    pass


class LOCATIONS:
    """These locations define where a reportable exception occurred and
    can be used to display it in the appropriate location in the report.
    """
    BLAST = 1.0
    BOLD_ID_ENGINE = 1.10
    BOLD_TAXA = 1.11
    DATABASE_COVERAGE = 5.0
    DATABASE_COVERAGE_NO_GBIF_RECORD = 5.01
    DATABASE_COVERAGE_NO_TAXID = 5.02
    DB_COVERAGE_TARGET = 5.1
    DB_COVERAGE_RELATED = 5.2
    DB_COVERAGE_RELATED_COUNTRY = 5.3


def write(location, msg, exc, query_dir=None, context=None):
    """Write a non-fatal error to file for later recall.

    location: display the error message in an appropriate
              location in the report.
    msg: provide the detailed information about the errors.
    exception: optional, can be None.
    query_dir: the location of error output file
               if the error happens in any specific query.
    data: optional, the context data of the error.
    """
    parent = query_dir or config.output_dir
    next_path = parent / config.ERRORS_DIR / 'next.txt'
    if next_path.exists():
        i = int(next_path.read_text())
    else:
        next_path.parent.mkdir(parents=True, exist_ok=True)
        i = 1
    next_path.write_text(str(i + 1))
    path = parent / config.ERRORS_DIR / f'{i}.json'
    with path.open('w') as f:
        json.dump({
            "location": location,
            "message": msg,
            "exception": str(exc),
            "context": context,
        }, f)


def read(query_dir=None, index=None):
    """Read all error files from given error directory."""
    errors = {}
    parent = query_dir or config.output_dir
    for path in (parent / config.ERRORS_DIR).glob('*.json'):
        with path.open() as f:
            errors[int(path.stem)] = json.load(f)
    if index:
        return errors.get(index)
    return errors


def report(path: Path, query_dir=None, location_min=0, location_max=999):
    """Write the specified errors to the given path for reporting."""
    errors = [
        x for x in read(query_dir=query_dir).values()
        if x['location'] >= location_min
        and x['location'] <= location_max
    ]
    with path.open('w') as f:
        for error in errors:
            f.write(
                f'Message: {error["message"]}\n'
                f'Exception: {error["exception"]}\n\n'
            )
