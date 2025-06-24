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

## BLAST - parsing the XML output

BLAST search is performed using a local (meaning run on the same machine as the workflow) BLASTN from the [NCBI BLAST+ toolkit](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html); the version is specified in the workflow report. This command-line BLASTN process produces a series of alignments for each query sequence, with each alignment relating to a BLAST "hit" against a sequence in the reference database.

<p class="alert alert-info">
    This step forks a single BLAST result into a series of query directories. From here each query's results are analysed in parallel.
</p>

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
            <td class="text-center p-3">
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
            <td class="text-center p-3">
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
            <td class="text-center p-3">
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
            <td class="text-center p-3">
                \(
                    \frac{\text{alignment length}}{\text{query length}}
                \)
            </td>
        </tr>
    </tbody>
</table>

## BLAST - Extracting taxonomic metadata

BLAST results do not include structured taxonomic information. This data is extracted for each BLAST hit subject using [taxonkit](https://bioinf.shenwei.me/taxonkit/), a command-line tool which can retrieve taxonomic records from NCBI's [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) archive. Taxids for each hit are extracted from the local BLAST database using `blastdbcmd`, another tool in the BLAST+ suite. This results in the following fields being collected for each hit:

- Taxid
- Domain
- Superkingdom
- Kingdom
- Phylum
- Class
- Order
- Family
- Genus
- Species

## BOLD - submitting sequences to ID Engine

When the workflow is run in `--bold` mode, the search method changes to use the BOLD ID Engine through the BOLD API (http://v4.boldsystems.org/index.php/Ids_xml). Since BOLD requires query DNA sequences that are correctly orientated (i.e. not antisense), we attempt to orientate the query sequences before submission. Query sequences are then submitted to the ID Engine API on-by-one. BOLD then returns a set of match statistics similar to BLAST for each query.

<p class="alert alert-info">
    This step forks each BOLD into a series of query directories. From here each query's results are analysed in parallel.
</p>

### Sequence orientation

Each DNA sequence is translated all three translation frames in both the forward and reverse directions. This results in six translated amino acid sequences for each query in frames `1`, `2`, `3`, `-1`, `-2`, `-3`.

To orientate each query sequence, we then use the `hmmsearch` tool (part of the [HMMER suite](http://eddylab.org/software/hmmer/Userguide.pdf)) locally to determine whether any the translation frames contain any of the following HMM profiles:

- `pf00115.hmm` - [Cytochrome C and Quinol oxidase polypeptide I](https://www.ebi.ac.uk/interpro/entry/pfam/PF00115/)
- `pf00116.hmm` - [Cytochrome C oxidase subunit II, periplasmic domain](https://www.ebi.ac.uk/interpro/entry/pfam/PF00116/)
- `pf02790.hmm` - [Cytochrome C oxidase subunit II, transmembrane domain](https://www.ebi.ac.uk/interpro/entry/pfam/PF02790/)

A match is accepted when the E-value is below `1e-5`. The first frame which is predicted to encode one of these domains dictates the orientation that will then be submitted to BOLD. For query sequences with no matches, both the forward and reverse orientations are submitted to BOLD and the one which returns hit(s) is assumed to be in the correct orientation (the other orientation's result is discarded).

### Submitting to ID Engine

Orientated query sequences are submitted to the ID Engine API sequentially, and the requests run in parallel to increase throughput.
The following data are parsed directly from the API response:

- Query title
- Query length
- Query frame
- Query sequence
- Hits:
    - Hit identifier (BOLD ID)
    - Hit sequence description
    - Hit taxonomic identification (species)
    - Hit similarity (used in place of identity)
    - Hit URL (a link to the record on [https://boldsystems.org](https://boldsystems.org))
    - Hit nucleotide sequence
    - Hit collectors (BOLD database submitter(s))

### Requesting additional metadata

For each hit subject, additional metadata are then requested from the "Full data retrieval" BOLD API endpoint:

- Accession (GenBank accession)
- Phylum
- Class
- Order
- Family
- Genus
- Species

The above fields are then used to fetch a kingdom classification (not included in BOLD response data) from the GBIF API.

## Assigning taxonomic identity



## Assessment of supporting publications


## Assessment of database coverage


## Report generation
