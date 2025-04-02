"""Provide an interface to the BOLD API for search and metadata retrieval."""

import csv
import logging
import requests
from Bio import SeqIO
from pathlib import Path
from xml.etree import ElementTree

from src.utils import errors
from src.utils.throttle import ENDPOINTS, Throttle
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)
ID_ENGINE_URL = "http://v4.boldsystems.org/index.php/Ids_xml"
FULL_DATA_URL = "http://v4.boldsystems.org/index.php/API_Public/combined"
SEQUENCE_DATA_URL = "http://v4.boldsystems.org/index.php/API_Public/sequence"
MIN_HTTP_CODE_ERROR = 400


class BoldSearch:
    """Fetch metadata for given taxa from the BOLD API."""
    def __init__(self, fasta_file: Path):
        self.fasta_file = fasta_file
        raw_hits = self._bold_sequence_search(fasta_file)
        if not any(raw_hits.values()):
            logger.info("No hits from BOLD ID Engine for FASTA query.")
        self.hits = self._fetch_hit_metadata(raw_hits)
        self.taxa = self._extract_taxa(self.hits)
        self.records = self._fetch_records(self.taxa)

    @property
    def taxon_count(self) -> dict[str, int]:
        """Return a count of records for each taxon."""
        taxa_counts = {}
        for record in self.records:
            taxon = record.get("species_name", "")
            taxa_counts[taxon] = taxa_counts.get(taxon, 0) + 1
        return taxa_counts

    @property
    def taxon_collectors(self) -> dict[str, list[str]]:
        """Return a dictionary of taxa and related collectors."""
        taxa_collectors = {}
        for record in self.records:
            taxon = record.get("species_name", "")
            collector = record.get("collectors", "")
            if taxon:
                if taxon not in taxa_collectors:
                    taxa_collectors[taxon] = set()
                if collector:
                    taxa_collectors[taxon].update(collector.split(","))
        return taxa_collectors

    @property
    def taxon_taxonomy(self) -> dict[str, dict[str, str]]:
        """Build a taxonomy dictionary."""
        taxonomy_dict = {}
        for record in self.records:
            species_name = record.get("species_name", "")
            if species_name:
                taxonomy_dict[species_name] = {
                    "phylum": record.get("phylum_name", ""),
                    "class": record.get("class_name", ""),
                    "order": record.get("order_name", ""),
                    "family": record.get("family_name", ""),
                    "genus": record.get("genus_name", ""),
                    "species": species_name,
                }
        return taxonomy_dict

    def _read_sequence_from_fasta(
        self,
        fasta_file: Path,
    ) -> list[SeqIO.SeqRecord]:
        """Read sequence from fasta file."""
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
        return sequences

    def _bold_sequence_search(
        self,
        fasta_file: Path,
        db: str = "COX1_SPECIES_PUBLIC",
    ) -> dict[str, list[dict[str, any]]]:
        """Submit a sequence search request to BOLD API with throttling."""
        hits = {}
        throttle = Throttle(ENDPOINTS.BOLD)
        sequences = self._read_sequence_from_fasta(fasta_file)
        logger.debug(
            f"Submitting {len(sequences)} sequences to BOLD ID Engine..."
        )
        for i, sequence in enumerate(sequences):
            logger.debug(
                f"Submitting sequence {i + 1}/{len(sequences)}: {sequence.id}"
            )
            params = {
                "sequence": str(sequence.seq),
                "db": db
            }
            response = throttle.with_retry(
                requests.get,
                args=[ID_ENGINE_URL],
                kwargs={"params": params})
            if response.status_code >= MIN_HTTP_CODE_ERROR:
                msg = (
                    f"Error for sequence {i + 1}: HTTP {response.status_code}"
                )
                logger.error(msg)
                errors.write(
                    errors.LOCATIONS.BOLD_ID_ENGINE,
                    msg,
                    None,
                    data={
                        "query_ix": i,
                        "response": response.text,
                    },
                )
                hits[sequence.id] = []
                continue

            root = ElementTree.fromstring(response.text)
            sequence_hits = []
            for match in root.findall("match"):
                result = {
                    "hit_id": match.find("ID").text,
                    "sequence_description": (
                        match.find("sequencedescription").text
                    ),
                    "database": match.find("database").text,
                    "citation": match.find("citation").text,
                    "taxonomic_identification": (
                        match.find("taxonomicidentification").text
                    ),
                    "similarity": match.find("similarity").text,
                    "specimen": match.find("specimen").text,
                    "url": match.find("specimen/url").text,
                    "country": (
                        match.find("specimen/collectionlocation/country").text
                    ),
                    "latitude": (
                        match.find(
                            "specimen/collectionlocation/coord/lat"
                        ).text
                        if match.find(
                            "specimen/collectionlocation/coord/lat"
                        ) is not None else ""
                    ),
                    "longitude": (
                        match.find(
                            "specimen/collectionlocation/coord/lon"
                        ).text
                        if match.find(
                            "specimen/collectionlocation/coord/lon"
                        ) is not None else ""
                    ),
                }
                sequence_hits.append(result)
            hits[sequence.id] = sequence_hits
        return hits

    def _fetch_hit_metadata(self, hits) -> dict[str, list]:
        """Fetch metadata by calling BOLD public API with accessions."""
        def _fetch_batch(batch_ids):
            params = {
                "ids": "|".join(batch_ids),
                "format": "tsv",
            }
            throttle = Throttle(ENDPOINTS.BOLD)
            response = throttle.with_retry(
                requests.get,
                args=[FULL_DATA_URL],
                kwargs={"params": params},
            )
            if response.status_code >= MIN_HTTP_CODE_ERROR:
                response.raise_for_status()
            lines = response.text.splitlines()
            reader = csv.DictReader(lines, delimiter="\t")
            return {
                row.pop("processid"): row
                for row in reader
            }

        hit_record_ids = [
            hit["hit_id"]
            for query_hits in hits.values()
            for hit in query_hits
        ]
        logger.info(
            f"Fetching sequences for {len(hit_record_ids)} hit records..."
        )
        metadata = {}
        batch_size = 25

        with ThreadPoolExecutor(max_workers=50) as executor:
            futures = []
            for i in range(0, len(hit_record_ids), batch_size):
                batch = hit_record_ids[i:i + batch_size]
                futures.append(executor.submit(_fetch_batch, batch))
            for future in as_completed(futures):
                metadata.update(future.result())

        hits_with_metadata = {
            query_title: [
                {**hit, **metadata.get(hit["hit_id"], {})}
                for hit in query_hits
            ]
            for query_title, query_hits in hits.items()
        }
        self.raw_hits = None

        return hits_with_metadata

    def _extract_taxa(self, hits: dict[str, list[dict]]) -> list[str]:
        """Extract taxa (taxonomic_identification)."""
        taxa = []
        # print(f"Hits Result in _extract_taxa: {hits}")
        for matches in hits.values():
            for hit in matches:
                taxonomic_identification = hit.get("taxonomic_identification")
                if (
                    taxonomic_identification
                    and taxonomic_identification not in taxa
                ):
                    taxa.append(taxonomic_identification)
        return taxa

    def _fetch_records(self, taxa: list[str]) -> list[dict]:
        """Fetch accessions by calling BOLD public API with Taxa
        and throttling, and save accessions to a .tsv file."""
        records = []
        if not taxa:
            logger.info("No taxa to fetch BOLD metadata.")
            return records
        taxa_param = "|".join(taxa)
        params = {
            "taxon": taxa_param,
            "format": "tsv"
        }
        throttle = Throttle(ENDPOINTS.BOLD)
        response = throttle.with_retry(
                requests.get,
                args=[FULL_DATA_URL],
                kwargs={"params": params})
        if response.status_code <= MIN_HTTP_CODE_ERROR:
            lines = response.text.splitlines()
            if not lines:  # Check if 'lines' is empty
                msg = "Empty response received from BOLD API"
                logger.error(msg)
                errors.write(
                    errors.LOCATIONS.BOLD_TAXA,
                    msg,
                    None,
                    data={"taxa": taxa},
                )
                return records

            headers = lines[0].split("\t")
            for line in lines[1:]:
                row = dict(zip(headers, line.split("\t")))
                records.append(row)
            logger.info(
                f"Records fetched successfully: {len(records)} records."
            )
        else:
            msg = (
                f"Error HTTP {response.status_code}: {response.text}"
            )
            logger.error(msg)
            errors.write(
                errors.LOCATIONS.BOLD_TAXA,
                msg,
                None,
                data={
                    "taxa": taxa,
                    "response": response.text,
                },
            )

        return records
