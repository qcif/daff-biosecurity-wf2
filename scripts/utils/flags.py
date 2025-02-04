"""Describe how flags are used to record outcomes."""

FLAG_DETAILS = {
    1: {
        "name": "Interspecies diversity",
        "explanation": {
            "A": "1 candidate species matched with high stringency"
                 " (identity >= 98.5%)",
            "B": "2-3 candidate species matched with high stringency"
                 " (identity >= 98.5%)",
            "C": ">3 candidate species matched with high stringency"
                 " (identity >= 98.5%)",
            "D": "At least one candidate species matched with moderate"
                 " stringency (identity >=93.5% and <98.5%)",
            "E": "No candidate species matched",
        },
    },
    2: {
        "name": "Taxa of interest",
        "explanation": {
            "A": "Taxon of interest detected",
            "B": "Taxon of interest NOT detected",
        },
    },
    4: {
        "name": "Intraspecies diversity",
        "explanation": {
            "A": "There are diverse results (>5 sources) for this species"
                 " in the results",
            "B": "There are limited results (1-5) for this species in the"
                 " results",
        },
    },
    5.1: {
        "name": "Database representation for target species",
        "explanation": {
            "A": "This locus for this species is well represented in reference"
                 " database (>5 entries) ",
            "B": "This locus for this species has limited representation in"
                 " reference database (<=5 entries) ",
            "C": "This locus for this species is not present in reference"
                 " database (0 entries) ",
        },
    },
    5.2: {
        "name": "Database representation for related species",
        "explanation": {
            "A": ">90% of species have sequence for this locus",
            "B": "10-90% of species have sequence for this locus",
            "C": "<=10% of species have sequence for this locus",
        },
    },
    5.3: {
        "name": "Database representation for related species in country of"
                " origin",
        "explanation": {
            "A": "All species in genus with observations in country of origin"
                 " have"
                 " sequence for this locus",
            "B": "Not all species in genus with observations in country of"
                 " origin have"
                 " sequence for this locus",
            "C": "No data: No species in genus with observations in country of"
                 " origin",
        },
    },
    6: {
        "name": "Intraspecies diversity",
        "explanation": {
            "A": "Candidate species grouped within single clade",
            "B": "The query sequence clusters as sister to a species-specific"
                 " clade or as a single branch in a tree",
            "C": "Genotype diversity grouped into multiple diverse groups",
        },
    },
    7: {
        "name": "Preliminary identification confirmation",
        "explanation": {
            "A": "Candidate species is consistent with preliminary morphology"
                 " ID",
            "B": "Candidate species is not consistent with preliminary"
                 " morphology ID",
        },
    },
}
