While the HTML reports from this workflow aim to be self-documenting, there are many analytical details which the user may wish to understand in further detail. The purpose of this document is to provide this enhanced understanding of what exactly happens during the analysis, and how outcomes are derived from the resulting data.

- To set up and run the workflow, visit the Nextflow workflow repository: [qcif/nf-daff-biosecurity-wf2](https://github.com/qcif/nf-daff-biosecurity-wf2).
- For analysis code, see the Python modules repository (used in the above workflow): [qcif/daff-biosecurity-wf2](https://github.com/qcif/daff-biosecurity-wf2)

## Reference data

### BLAST

The reference data used by the workflow depends entirely on the deployment - ask your platform administrator if you are unsure.
For the BLAST version of the workflow, the reference data will be a BLAST database of sequence records that is held on the analysis server - by default this is NCBI's "Core Nt" database, but in future it may be possible to choose a different reference database. The HTML report should specify the database name in the database coverage report (admins can set this manually with `BLAST_DATABASE_NAME`).

### BOLD

By default this is set to `COX1_SPECIES_PUBLIC` (admins can override this by setting `BOLD_DATABASE`).

## Input files

### FASTA file

A FASTA file containing sample sequences to be analysed. Multiple sequences per sample can be used, but the FASTA header for each sequence must be unique and match an entry in the `metadata.csv` input. The following constraints apply to this input:

- Seq IDs must be unique
- Seq IDs must match `metadata.csv` input
- Max 150 sequences
- Minimum seq length: `20nt`
- Max length of any sequence: `3000nt`
- All residues must be valid nucleotide (ambiguous IUPAC DNA: `ATGCRYSWKMBDHVN`)

### Metadata CSV file

This file provides metadata for each query sequence, with the following fields:

| Field             | Required | Description                                                                                                         |
|-------------------|----------|---------------------------------------------------------------------------------------------------------------------|
| sample_id         | yes      | Must match the header of one FASTA sequence                                                                         |
| locus             | yes      | Must be in the [list of allowed Loci](./allowed-loci.html) or `NA` for virus or BOLD queries                                               |
| preliminary_id    | yes      | A suggested taxonomic identity based on sample morphology                                                           |
| taxa_of_interest  | no       | A pipe-delimited list of taxa to be evaluated against the sample. Can be at rank species, genus, family, order, class, phylum, kingdom or domain. |
| country           | no       | The sample country of origin                                                                                        |
| host              | no       | The host or commodity that the sample was extracted from                                                            |

<p class="alert alert-info">
  From <code>v1.1</code> arbitrary fields in this file will be displayed in the workflow report
</p>

## BLAST - parsing the BLAST result

BLAST search is performed using a local (meaning run on the same machine as the workflow) BLASTN from the [NCBI BLAST+ toolkit](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html); the version is specified in the workflow report. This command-line BLASTN process produces a series of alignments for each query sequence, with each alignment relating to a BLAST "hit" against a sequence in the reference database.

### Extracted values

The following values are extracted verbatim from BLAST XML fields:

- Hit identifier (GenBank ID)
- Hit definition (GenBank record title)
- Hit NCBI accession
- Hit subject length (nt)
- High-scoring pairs (HSPs; each represents a segment of the alignment):
    - bitscore
    - evalue
    - identity
    - query_start
    - query_end
    - subject_start
    - subject_end
    - alignment_length

### Calculated values

These values are not present in the BLAST XML and are calculated from the extracted values:

<table class="table table-striped">
    <thead>
        <th>Value</th>
        <th>Description</th>
        <th>Equation</th>
    </thead>
    <tbody>
        <tr>
            <td>
                Alignment length
            </td>
            <td>
                The total non-overlapping length of all HSPs.
            </td>
            <td>
            </td>
        </tr>
        <tr>
            <td>
                Hit bitscore
            </td>
            <td>
                A score which takes into account both alignment strength and length. Calculated as the sum of bitscores across all HSPs.
            </td>
            <td>
                \( \sum_{HSP \in \text{HSPs}} \text{bits}(HSP) \)
            </td>
        </tr>
        <tr>
            <td>
                Hit E-value
            </td>
            <td>
                An expression of probability that the alignment occurred due to random chance, often expressed as an exponent to distinguish between very low numbers. If there is only one HSP, the `hsp.evalue` will be used. Otherwise, a formula is used.
            </td>
            <td>
                \( \text{ess} \cdot 2^{-\sum_{HSP \in \text{hit.HSPs}} \text{bits}(HSP)} \)
                <br>
                <em>Where <code>ess</code> is the effective search space specified in the BLAST XML output.</em>
            </td>
        </tr>
        <tr>
            <td>
                Hit identity
            </td>
            <td>
                The proportion of nucleotides which match between query and subject in the alignment. This is calculated as the weighted identity of HSPs (high-scoring pairs), clipped to a maximum of 1.
            </td>
            <td>
                \(
                    \frac{\sum_{HSP \in \text{HSPs}} \text{identities}(HSP)}{\sum_{HSP \in \text{HSPs}} \text{alignment length}(HSP)}
                \)
            </td>
        </tr>
        <tr>
            <td>
                Query coverage
            </td>
            <td>
                The proportion of the query sequence that is covered by the alignment with the reference sequence.
            </td>
            <td>
                \(
                    \frac{\text{alignment length}}{\text{query length}}
                \)
            </td>
        </tr>
    </tbody>
</table>


## BOLD - submitting sequences to ID Engine


## BLAST - Extracting taxonomic metadata


## Assigning taxonomic identity


## Assessment of supporting publications


## Assessment of database coverage


## Report generation
