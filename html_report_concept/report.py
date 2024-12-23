"""Create an HTML report from template."""

import os
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

PWD = Path(__file__).parent
TEMPLATE_DIR = PWD / 'templates'
OUTPUT_DIR = PWD / 'output'
j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))

# Render the template
template = j2.get_template('report.html')
rendered_html = template.render(title="Workflow report")

# Save the rendered HTML to a file
output_file = os.path.join(OUTPUT_DIR, 'report.html')
with open(output_file, 'w') as file:
    file.write(rendered_html)

print(f"HTML document has been created: {output_file}")
