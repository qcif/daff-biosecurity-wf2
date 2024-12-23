import base64
import os
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

# Paths
PWD = Path(__file__).parent
TEMPLATE_DIR = PWD / 'templates'
TEMPLATE_NAME = "bam-viewer.html"
OUTPUT_DIR = PWD / 'output'
BAM_FILE_PATH = PWD / "data/long-reads.bam"
BAI_FILE_PATH = PWD / "data/long-reads.bam.bai"
OUTPUT_HTML_PATH = OUTPUT_DIR / TEMPLATE_NAME


def encode_file_to_base64(file_path):
    with open(file_path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


# Check if BAM and BAI files exist
if not os.path.exists(BAM_FILE_PATH) or not os.path.exists(BAI_FILE_PATH):
    raise FileNotFoundError("Both BAM and BAI files are required.")

bam_base64 = encode_file_to_base64(BAM_FILE_PATH)
bai_base64 = encode_file_to_base64(BAI_FILE_PATH)

j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
template = j2.get_template(TEMPLATE_NAME)

rendered_html = template.render(
    bam_data=bam_base64,
    bai_data=bai_base64,
)

# Write output HTML
with open(OUTPUT_HTML_PATH, "w") as f:
    f.write(rendered_html)

print(f"BAM Viewer HTML generated: {OUTPUT_HTML_PATH}")
