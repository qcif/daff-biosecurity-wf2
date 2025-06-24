While the workflow report aims to be self-documenting, there are many analytical details which the user may wish to understand in further detail. The purpose of this document is to provide this enhanced understanding of what exactly happens during the analysis, and how outcomes are derived from the resulting data.

- To set up and run the workflow, visit the Nextflow workflow repository: [qcif/nf-daff-biosecurity-wf2](https://github.com/qcif/nf-daff-biosecurity-wf2).
- For analysis code, see the Python modules repository (used in the above workflow): [qcif/daff-biosecurity-wf2](https://github.com/qcif/daff-biosecurity-wf2)

## Interpretting the workflow report

(This should probably be in the workflow report!)

## Reference data

### BLAST

The reference data used by the workflow depends entirely on the deployment - ask your platform administrator if you are unsure.
For the BLAST version of the workflow, the reference data will be a BLAST database of sequence records that is held on the analysis server - by default this is `{{ config.REPORT.DATABASE_NAME }}`, but it is possible to run with a different reference database. The workflow report specifies the database name in the database coverage report (admins can set this manually with `BLAST_DATABASE_NAME`).

### BOLD

By default this is set to `{{ config.BOLD_DATABASE }}` (admins can override this by setting `BOLD_DATABASE`).

## Input files

### FASTA file

A FASTA file containing sample sequences to be analysed. Multiple sequences per sample can be used, but the FASTA header for each sequence must be unique and match an entry in the `metadata.csv` input. The following constraints apply to this input:

- Seq IDs must be unique
- Seq IDs must match `metadata.csv` input
- Maximum query sequences: `{{ config.INPUTS.FASTA_MAX_SEQUENCES }}`
- Minimum seq length: `{{ config.INPUTS.FASTA_MIN_LENGTH_NT }}nt`
- Max length of any sequence: `{{ config.INPUTS.FASTA_MAX_LENGTH_NT }}nt`
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
  From <code>v1.1</code> any additional fields in this file will be displayed in the workflow report
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

When the workflow is run in `--bold` mode, the search method changes to use the BOLD ID Engine through the BOLD API ([https://v4.boldsystems.org/resources/api](https://v4.boldsystems.org/resources/api)). Since BOLD requires query DNA sequences that are correctly orientated (i.e. not antisense), we attempt to orientate the query sequences before submission. Query sequences are then submitted to the ID Engine API on-by-one. BOLD then returns a set of match statistics similar to BLAST for each query.

<p class="alert alert-info">
    This step forks each BOLD into a series of query directories. From here each query's results are analysed in parallel.
</p>

### Sequence orientation

Each DNA sequence is translated all three translation frames in both the forward and reverse directions. This results in six translated amino acid sequences for each query in frames `1`, `2`, `3`, `-1`, `-2`, `-3`.

To orientate each query sequence, we then use the `hmmsearch` tool (part of the [HMMER suite](http://eddylab.org/software/hmmer/Userguide.pdf)) locally to determine whether any the translation frames contain any of the following HMM profiles:

- `pf00115.hmm` - [Cytochrome C and Quinol oxidase polypeptide I](https://www.ebi.ac.uk/interpro/entry/pfam/PF00115/)
- `pf00116.hmm` - [Cytochrome C oxidase subunit II, periplasmic domain](https://www.ebi.ac.uk/interpro/entry/pfam/PF00116/)
- `pf02790.hmm` - [Cytochrome C oxidase subunit II, transmembrane domain](https://www.ebi.ac.uk/interpro/entry/pfam/PF02790/)

A match is accepted when the E-value is below `{{ config.HMMSEARCH_MIN_EVALUE }}`. The first frame which is predicted to encode one of these domains dictates the orientation that will then be submitted to BOLD. For query sequences with no matches, both the forward and reverse orientations are submitted to BOLD and the one which returns hit(s) is assumed to be in the correct orientation (the other orientation's result is discarded).

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
The phylum given above is used to fetch a the associated phylum record from the [GBIF species search API](https://techdocs.gbif.org/en/openapi/v1/species#/Searching%20names/searchNames), and the kingdom field is extracted from the response data.

## Assigning taxonomic identity

This is a critical stage in the analysis. Hits returned from BLAST/BOLD are filtered and a list of candidate species is extracted from those hits. Filtering is applied as follows. For BOLD search, the process is identical, with similarity being used in place of identity.

- All hits which are below `{{ config.CRITERIA.ALIGNMENT_MIN_NT }}nt` AND `{{ config.CRITERIA.ALIGNMENT_MIN_Q_COVERAGE * 100 }}%` query coverage are excluded from the entire analysis. These are referred to as "filtered hits".
- The identity threshold for candidate hits is either:
    - `{{ config.CRITERIA.ALIGNMENT_MIN_IDENTITY_STRICT * 100 }}%` (if any filtered hits meet this threshold) - defined as a **STRONG MATCH**
    - OR `{{ config.CRITERIA.ALIGNMENT_MIN_IDENTITY * 100 }}%` - defined as a **MODERATE MATCH**
- Resulting candidate hits are then aggregated by species and the top hit per-species is identified. These species are what you see reported in the "Candidates" section of the workflow report.
- The "No. hits" shown for each candidate includes all filtered hits, not just the candidate hits.
- The "Median identity" shown in the "Candidate species" table is derived from the identity (%) of all filtered hits. If there is a wide distribution of hit identities, the median will be reduced. If the median drops below the candidate threshold, the badge will turn from green to yellow. If it drops to less than `{{ config.CRITERIA.MEDIAN_IDENTITY_WARNING_FACTOR * 100 }}%` of the threshold (i.e. `<{{ (config.CRITERIA.ALIGNMENT_MIN_IDENTITY_STRICT * 100 * config.CRITERIA.MEDIAN_IDENTITY_WARNING_FACTOR) | round(1) }}%` for a threshold of `{{ config.CRITERIA.ALIGNMENT_MIN_IDENTITY_STRICT * 100 }}%`), the badge will turn red.

The candidate species identified above are then cross-checked against the Preliminary Morphology ID and Taxa of Interest (if provided).

### Checking preliminary morphology ID

The PMI provided by the user is checked against each taxonomic rank from each candidate species. If the provided name matches any of the fields, this is regarded as a match. The user should be aware that in some edge-cases this can result in a mis-match due to the existence of taxonomic homonyms (*Morus*, for example, is both a plant and bird genus).

### Checking taxa of interest

This process is identical to that described above, except that a little more information is collected for display in the "Taxa of interest" section of the workflow report. For each taxon of interest, the best-scoring species that matches the taxonomy is reported. It is possible for a TOI to match multiple candidates, but only the top candidate will be shown.


## Phylogenetic analysis

Subject sequences are selected from filtered hits [extracted previously](#assigning-taxonomic-identity).
The selection process is a little complex, as it aims to strike a balance between reasonable coverage of the genetic diversity present in BLAST/BOLD hit subjects, while also trying to minimize the number of sequences that need to go through alignment and analysis. Building trees with 100+ sequences is SLOW and the resulting tree is often ugly, so we do our best to avoid that.

- Hits are collected in order of descending identity until both:
    - The identity of the hit drops to `{{ config.CRITERIA.PHYLOGENY_MIN_HIT_IDENTITY * 100 }}%`
    - AND at least {{ config.CRITERIA.PHYLOGENY_MIN_HIT_SEQUENCES }} hits have been collected

This means that candidate hits are always sampled, and some filtered hits are usually included too.

Next, if there are more than {{ config.CRITERIA.PHYLOGENY_MAX_HITS_PER_SPECIES }} hits for the species, a representative sample of hits for each species is taken:

- Filtered hits for the species are ordered by identity
- A systematic sample of {{ config.CRITERIA.PHYLOGENY_MAX_HITS_PER_SPECIES }} hits is taken, including the first and last hit

The aim is that a species with 200 hits of alignment identities between 90-95% would result in a sample of hits with identities similar to the following (assuming a bimodal distribution representing two taxonomic clades, and sample size of n=6):

- 90.0%
- 90.3%
- 90.5%
- 94.8%
- 94.8%
- 95.0%

<p class="alert alert-info">
    The workflow has a default value of <code>PHYLOGENY_MAX_HITS_PER_SPECIES=1000</code> so that the sampling described above is not implemented unless the workflow administrator decides to limit this number to constrain tree size. Setting this sample size too low would result in poor quality trees that may give a false impression of genetic diversity to the user.
</p>

A FASTA sequence is then written by extracting the nucleotide sequence from each of the selected hits, and adding the query sequence.
This FASTA file is then alignment with [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html), and then a tree is computed from the alignment with [FastME](http://www.atgc-montpellier.fr/fastme/).

## Assessment of supporting publications

An important drawback of searching against the large Non-redundant (Nr) BLAST database is that this database contains many sequences which are not very reputable. Since anyone can submit sequences to GenBank there are many sequences with incorrect taxonomic annotation. This analysis presents a measure of confidence in the integrity of the reference sequences supporting the conclusions. Are the candidate reference sequences supported by numerous publications? Great, that means that the taxonomic annotation has been corroborated by multiple studies. Do we have 5 reference sequences that were all submitted to genbank by the same author(s)? That casts some doubt over the integrity of the taxonomic annotation.

<p class="alert alert-warning">
    Even when the workflow report is "green" in all other sections of the report, caution is advised if there is only one independent publication corroborating the reference sequences. Further investigation may be required to confirm the credibility of the reference sequence source.
</p>

The analysis involves clusting of GenBank publication records based on the provided metadata (author, title, journal). An independent analysis is carried out per-candidate species:

1. For each [candidate hit](#assigning-taxonomic-identity), a list of publications is extracted from the corresponding GenBank record using the Entrez `efetch` API.
1. Each publication is annotated with:
    1. The NCBI accession number
    1. A `source tag`, which is a token derived (by stripping non-alphanumeric characters and converting to lowercase) from one of the following fields (the field with a value):
        1. Author list
        2. Publication title ("Direct submission" titles are ignored here as they are computed-generated and not a real publication)
        3. Journal
    1. Whether the record appears to have been automatically generated - determined by the presence of the string `##Genome-Annotation-Data-START##` in the GenBank record text.
1. Publications are then clustered into groups of "independent publications":
    1. If a publication `source tag` matches one in an existing group, it is added to that group
    1. Otherwise, it is added to a new group
    1. Records with no publications are allocated to a single group
1. The groups are shown as "independent sources" in the "Publications supporting reference sequences" section of the workflow report.

## Assessment of database coverage



## Report generation
