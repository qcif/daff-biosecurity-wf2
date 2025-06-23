"""Render documentation in Markdown format to HTML."""

from pathlib import Path

import markdown2
from jinja2 import Environment, FileSystemLoader

ROOT_DIR = Path(__file__).resolve().parents[2]
DOCS_SRC_DIR = ROOT_DIR / 'docs/src'
DOCS_DEST_DIR = ROOT_DIR / 'docs'


def main():
    DOCS_DEST_DIR.mkdir(parents=True, exist_ok=True)
    for md_file in DOCS_SRC_DIR.glob('*.md'):
        dest_filename = md_file.stem + '.html'
        dest_path = DOCS_DEST_DIR / dest_filename
        md_content = md_file.read_text(encoding='utf-8')
        html_body = markdown2.markdown(md_content, extras={
            "tables": True,
            "code-friendly": True,
            "html-classes": {
                'table': 'table table-striped',
            },
        })
        j2 = Environment(
            loader=FileSystemLoader(str(ROOT_DIR / 'docs/src')),
            autoescape=True,
        )
        html_doc = j2.get_template('header.html').render(
            body=html_body,
        )
        dest_path.write_text(html_doc, encoding='utf-8')
        print(f"Rendered {md_file} -> {dest_path}")


if __name__ == '__main__':
    main()
