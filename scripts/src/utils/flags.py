"""Describe how flags are used to record outcomes."""

import logging

from .config import Config

logger = logging.getLogger(__name__)
config = Config()


class Flag:
    """A flag to record discrete analytical outcomes.

    If a target is provided, it should be the target species name, with
    target_type being one of ['candidate', 'pmi', 'toi'].
    """

    def __init__(
        self,
        flag_id,
        value: str = None,
        target: str = None,
        target_type: str = None,
    ):
        self.flag_id = flag_id
        self.value = value
        self.target = target
        self.target_type = target_type

    def __str__(self):
        return f"{self.flag_id}{self.value}"

    def __repr__(self):
        return f"{self.flagCRITERIA.ALIGNMENT_MIN_id}{self.value}"

    def to_json(self):
        data = {
            "flag_id": self.flag_id,
            "value": self.value,
            "target": self.target,
            "level": self.level,
            "outcome": self.outcome,
            "explanation": self.explanation,
            "bs-class": self.bs_class,
        }

        return data

    @property
    def name(self):
        return FLAG_DETAILS[self.flag_id]["name"]

    @property
    def explanation(self):
        return FLAG_DETAILS[self.flag_id]["explanation"][self.value]

    @property
    def outcome(self):
        return FLAG_DETAILS[self.flag_id]["outcome"][self.value]

    @property
    def level(self):
        """Return the warning level for the given value."""
        return FLAG_DETAILS[self.flag_id]['level'][self.value]

    @property
    def bs_class(self):
        """Return the bootstrap css class for self.level."""
        level = self.level
        if level == 0:
            return "secondary"
        if level == 1:
            return "success"
        if level == 2:
            return "warning"
        return "danger"

    @classmethod
    def read(cls, query, as_json=False):
        """Read flags from CSV file."""
        pattern = (
            config.get_query_dir(query)
            / config.FLAG_FILE_TEMPLATE.format(identifier="*")
        )

        flags = {}
        for path in pattern.parent.glob(pattern.name):
            if path.is_file():
                flag_id, target, target_type = cls._parse_flag_filename(path)
                value = path.read_text().strip()
                if not as_json:
                    value = Flag(flag_id, value=value, target=target)
                if target:
                    flags[flag_id] = flags.get(flag_id, {})
                    if target_type:
                        flags[flag_id][target_type] = flags[flag_id].get(
                            target_type, {})
                        flags[flag_id][target_type][target] = value
                    else:
                        flags[flag_id][target] = value
                else:
                    flags[flag_id] = value
        for flag, data in flags.items():
            if flag.startswith('5'):
                for ttype in [TARGETS.CANDIDATE, TARGETS.PMI, TARGETS.TOI]:
                    if ttype not in data:
                        data[ttype] = {}
        return flags

    def _parse_flag_filename(path):
        """Parse flag filename to extract flag_id, target, target_type."""
        stem = path.stem
        target_type = None
        if '-' in stem:
            flag_id, target = stem.split("-", 1)
            if '-' in target:
                target_type, target = target.split("-")
            target = target.replace("_", " ")
        else:
            flag_id = stem
            target = None
        return flag_id, target, target_type

    @classmethod
    def write(
        cls,
        query_dir,
        flag_id,
        value,
        target=None,
        target_type=None,
    ):
        """Write flags value to JSON file.

        target: the target taxon name
        target_type: one of TARGETS.[candidate, pmi, toi]
        """
        identifier = flag_id
        target_str = ''
        if target:
            type_str = f"{target_type}-" if target_type else ''
            target_str = f"-{type_str}{target}".replace(' ', '_')
            identifier += target_str
            target_str = ' ' + target_str.strip('-')
        path = query_dir / config.FLAG_FILE_TEMPLATE.format(
            identifier=identifier)
        with path.open('w') as f:
            f.write(value)
        logger.info(f"Flag {flag_id}{value}{target_str} written to {path}")


class TARGETS:
    CANDIDATE = "candidate"
    PMI = "pmi"
    TOI = "toi"


class FLAGS:
    """Flags for reporting outcomes."""
    POSITIVE_ID = '1'
    TOI = '2'
    SOURCES = '4'
    DB_COVERAGE_TARGET = '5.1'
    DB_COVERAGE_RELATED = '5.2'
    DB_COVERAGE_RELATED_COUNTRY = '5.3'
    INTRASPECIES_DIVERSITY = '6'
    PMI = '7'
    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"


FLAG_DETAILS = config.read_flag_details_csv()
