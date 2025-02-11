"""Describe how flags are used to record outcomes."""

import json
import logging

from utils.config import Config

logger = logging.getLogger(__name__)

config = Config()


class FLAGS:
    """Flags for reporting outcomes."""
    POSITIVE_ID = '1'
    TOI = '2'
    SOURCES = '3'
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


class Flag:
    """A flag to record discrete analytical outcomes.

    If targets are provided, they should be a dictionary of TARGET.NAME:
    {ix: value}, where ix matches the order of the <TARGET.NAME> taxa.
    """

    def __init__(
        self,
        flag_id,
        value: str = None,
        targets: dict[str, dict[int, str]] = None,
    ):
        self.flag_id = flag_id
        self.value = value
        self.targets = targets

    def __str__(self):
        return f"{self.flag_id}{self.value}"

    def __repr__(self):
        return f"{self.flag_id}{self.value}"

    def to_json(self):
        data = {
            "flag_id": self.flag_id,
        }
        if self.targets:
            targets = {
                target: {
                    ix: {
                        "value": value,
                        "level": self.get_level(value),
                        "explanation": self.explanation(value),
                    }
                    for ix, value in value.items()
                }
                for target, value in self.targets.items()
            }
            data["targets"] = targets
        else:
            data["value"] = self.value
            data["level"] = self.get_level()
            data["explanation"] = self.explanation()
        return data

    @property
    def name(self):
        return FLAG_DETAILS[self.flag_id]["name"]

    def explanation(self, value=None):
        return FLAG_DETAILS[self.flag_id]["explanation"][value or self.value]

    def value_for_target(self, target, index=None, max_only=False):
        if target == TARGETS.PMI:
            return self.targets.get(target)
        if index is not None:
            try:
                return self.targets.get(target)[index]
            except KeyError:
                raise ValueError(f"Index {index} is out of range for flag"
                                 f" {self.flag_id} target {target}.")
        if max_only:
            return sorted(self.targets.get(target).values())[-1]
        raise ValueError(f"You must provide either index or max_value when"
                         f" requesting flag for target {target}.")

    def get_level(self, value=None):
        """Return the warning level for the given value."""
        return FLAG_DETAILS[self.flag_id]['level'][value or self.value]

    @staticmethod
    def read_flags(query_ix, as_json=False):
        """Read flags from CSV file."""
        path = config.get_query_dir(query_ix) / config.FLAGS_JSON
        if not path.exists():
            return {}
        with path.open() as f:
            if as_json:
                return json.load(f)
            return {
                flag_id: (
                    Flag(flag_id, value=data)
                    if isinstance(data, str)
                    else
                    Flag(flag_id, targets=data['targets'])
                )
                for flag_id, data in json.load(f).items()
            }

    @classmethod
    def write(
        cls,
        query_ix,
        flag_id,
        value,
        target=None,
        target_ix=None,
    ):
        """Write flags value to JSON file."""
        path = config.get_query_dir(query_ix) / config.FLAGS_JSON
        flags = cls.read_flags(query_ix, as_json=True)
        if target:
            flag = flags.get(flag_id, {})
            targets_dict = flag.get('targets', {})
            target_dict = targets_dict.get(target, {})
            target_dict[target_ix] = value
            targets_dict[target] = target_dict
            flag['targets'] = targets_dict
            flags[flag_id] = flag
        else:
            flags[flag_id] = value
        with path.open('w') as f:
            json.dump(flags, f, indent=2)
        logger.info(f"Flag {flag_id}{value} written to {config.FLAGS_CSV}")


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
    },
    FLAGS.TOI: {
        "name": "Taxa of interest",
        "explanation": {
            FLAGS.A: "Taxon of interest detected",
            FLAGS.B: "Taxon of interest NOT detected",
        },
        "level": {
            FLAGS.A: 1,
            FLAGS.B: 3,
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
