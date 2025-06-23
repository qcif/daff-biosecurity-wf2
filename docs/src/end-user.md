# User documentation

<p class="lead">DAFF Taxonomic assignment workflow</p>

While the HTML reports from this workflow aim to be self-documenting, there are many analytical details which the user may wish to understand in further detail. The purpose of this document is to provide this enhanced understanding of what exactly happens during the analysis, and how outcomes are derived from the resulting data.

To set up and run the workflow, visit the Nextflow workflow repository: [qcif/nf-daff-biosecurity-wf2](https://github.com/qcif/nf-daff-biosecurity-wf2).

For analysis code, see the Python modules repository (used in the above workflow): [qcif/daff-biosecurity-wf2](https://github.com/qcif/daff-biosecurity-wf2)

## Reference data

### BLAST

The reference data used by the workflow depends entirely on the deployment - ask your platform administrator if you are unsure.
For the BLAST version of the workflow, the reference data will be a BLAST database of sequence records that is held on the analysis server - by default this is NCBI's "Core Nt" database, but in future it may be possible to choose a different reference database. The HTML report should specify the database name in the database coverage report (admins can set this manually with `BLAST_DATABASE_NAME`).

### BOLD

By default this is set to `COX1_SPECIES_PUBLIC` (admins can override this by setting `BOLD_DATABASE`).

## Validation of user inputs


## BLAST - parsing the BLAST result


## BOLD - submitting sequences to ID Engine


## BLAST - Extracting taxonomic metadata


## Assigning taxonomic identity


## Assessment of supporting publications


## Assessment of database coverage


## Report generation
