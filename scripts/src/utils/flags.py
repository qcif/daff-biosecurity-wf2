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
        return f"{self.flag_id}{self.value}"

    def to_json(self):
        data = {
            "flag_id": self.flag_id,
            "value": self.value,
            "target": self.target,
            "level": self.get_level(),
            "outcome": self.outcome(),
            "explanation": self.explanation(),
            "bs-class": self.get_bs_class(),
        }

        return data

    @property
    def name(self):
        return FLAG_DETAILS[self.flag_id]["name"]

    def explanation(self):
        return FLAG_DETAILS[self.flag_id]["explanation"][self.value]

    def outcome(self):
        return FLAG_DETAILS[self.flag_id]["outcome"][self.value]

    def get_level(self):
        """Return the warning level for the given value."""
        return FLAG_DETAILS[self.flag_id]['level'][self.value]

    def get_bs_class(self):
        """Return the warning level for the given value."""
        level = self.get_level(self.value)
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
                    flags[flag_id][target_type] = flags[flag_id].get(
                        target_type, {})
                    flags[flag_id][target_type][target] = value
                else:
                    flags[flag_id] = value
        return flags

    def _parse_flag_filename(path):
        """Parse flag filename to extract flag_id, target, target_type."""
        stem = path.stem.split('flag_', 1)[1]
        if '[' in stem:
            flag_id, target = stem.split("[")
            target_type, target = target.split("-")
            target = target.replace("]", "").replace("_", " ")
        else:
            flag_id = stem
            target = None
            target_type = None
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
        target_type: one of ['candidate', 'pmi', 'toi']
        """
        identifier = (
            f"{flag_id}[{target_type}-{target}]"
            if target
            else flag_id
        ).replace(' ', '_')
        path = query_dir / config.FLAG_FILE_TEMPLATE.format(
            identifier=identifier)
        with path.open('w') as f:
            f.write(value)
        target_str = f" (target {target_type}:{target})" if target else ""
        logger.info(f"Flag {flag_id}{value}{target_str} written to {path}")


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


class TARGETS:
    CANDIDATE = "candidate"
    PMI = "pmi"
    TOI = "toi"


FLAG_DETAILS = {
    FLAGS.POSITIVE_ID: {
        "name": "Interspecies diversity",
        "explanation": {
            FLAGS.A: "1 candidate species matched with high stringency"
                " (identity >= 98.5%)",
            FLAGS.B: "2-3 candidate species matched with high stringency"
                " (identity >= 98.5%)",
            FLAGS.C: ">3 candidate species matched with high stringency"
                " (identity >= 98.5%)",
            FLAGS.D: "At least one candidate species matched with moderate"
                " stringency (identity >=93.5% and <98.5%)",
            FLAGS.E: "No candidate species matched",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
            FLAGS.C: 3,
            FLAGS.D: 2,
            FLAGS.E: 3,
        },
        "outcome": {
            FLAGS.A: "Positive species identification",
            FLAGS.B: "The analyst should attempt subjective species"
                     " identification at the genus level",
            FLAGS.C: "The analyst should attempt subjective species"
                     " identification at the genus level",
            FLAGS.D: "The analyst should attempt subjective species"
                     " identification at the genus level",
            FLAGS.E: "Identification not possible (potential unknown species)",
        }
    },
    FLAGS.TOI: {
        "name": "Taxa of interest",
        "explanation": {
            FLAGS.A: "Taxon of interest NOT detected",
            FLAGS.B: "Taxon of interest detected",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 3,
        },
        "outcome": {
            FLAGS.A: "None of the taxa of interest matched the"
                     " candidates species",
            FLAGS.B: "At least one of the taxa of interest matched the"
                     " candidates species",
        },
    },
    FLAGS.SOURCES: {
        "name": "Reference sequence source diversity",
        "explanation": {
            FLAGS.A: "There are diverse results (>5 sources) for this species"
            " in the results",
            FLAGS.B: "There are limited results (1-5) for this species in the"
            " results",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
        },
    },
    FLAGS.DB_COVERAGE_TARGET: {
        "name": "Database representation for target species",
        "explanation": {
            FLAGS.A: "This locus for this species is well represented in"
                " reference database (>5 entries) ",
            FLAGS.B: "This locus for this species has limited representation"
                " in reference database (<=5 entries) ",
            FLAGS.C: "This locus for this species is not present in reference"
                " database (0 entries) ",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
            FLAGS.C: 3,
        },
    },
    FLAGS.DB_COVERAGE_RELATED: {
        "name": "Database representation for related species",
        "explanation": {
            FLAGS.A: ">90% of species have sequence for this locus",
            FLAGS.B: "10-90% of species have sequence for this locus",
            FLAGS.C: "<=10% of species have sequence for this locus",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
            FLAGS.C: 3,
        },
    },
    FLAGS.DB_COVERAGE_RELATED_COUNTRY: {
        "name": "Database representation for related species in country of"
                " origin",
        "explanation": {
            FLAGS.A: "All species in genus with observations in country of"
                " origin have sequence for this locus",
            FLAGS.B: "Not all species in genus with observations in country of"
                " origin have"
                " sequence for this locus",
            FLAGS.C: "No data: No species in genus with observations in"
                " country of origin",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
            FLAGS.C: 0,
        },
    },
    FLAGS.INTRASPECIES_DIVERSITY: {
        "name": "Intraspecies diversity",
        "explanation": {
            FLAGS.A: "Candidate species grouped within single clade",
            FLAGS.B: "The query sequence clusters as sister to a"
                " species-specific clade or as a single branch in a tree",
            FLAGS.C: "Genotype diversity grouped into multiple diverse groups",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 2,
            FLAGS.C: 3,
        },
    },
    FLAGS.PMI: {
        "name": "Preliminary identification confirmation",
        "explanation": {
            FLAGS.A: "Candidate species is consistent with preliminary"
                " morphology ID",
            FLAGS.B: "Candidate species is not consistent with preliminary"
                " morphology ID",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 3,
        },
    },
}
