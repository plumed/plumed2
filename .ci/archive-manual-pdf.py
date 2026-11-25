#!/usr/bin/env python3
"""Package an unsupported manual's PDF for offline download on GitHub Pages.

Run on a prepared publication tree, or on a legacy doc-v2.N gh-pages checkout.
GitHub Pages cannot set X-Robots-Tag for PDFs. Keep the offline copy in a ZIP
instead of serving an indexable PDF, and update local HTML download links.
"""

import argparse
import posixpath
import re
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit
from zipfile import ZIP_DEFLATED, ZipFile


def archive_manual(root):
    pdf = root / "user-doc/manual.pdf"
    if not pdf.is_file():
        return

    # Finish and verify the download before removing the original PDF.
    archive = pdf.with_suffix(".pdf.zip")
    with ZipFile(archive, "w", compression=ZIP_DEFLATED) as output:
        output.write(pdf, arcname=pdf.name)
    with ZipFile(archive) as output:
        if output.read(pdf.name) != pdf.read_bytes():
            raise RuntimeError("PDF archive verification failed")

    for page in root.rglob("*.html"):
        parent = page.relative_to(root).parent.as_posix()

        def update_link(match):
            url = urlsplit(match[2])
            if url.scheme or url.netloc or url.path.startswith("/"):
                return match[0]
            target = posixpath.normpath(posixpath.join(parent, url.path))
            if target != "user-doc/manual.pdf":
                return match[0]
            # A PDF page fragment has no meaning on a ZIP download.
            href = urlunsplit(("", "", url.path + ".zip", url.query, ""))
            return match[1] + href + match[3]

        original = page.read_text(encoding="utf-8")
        updated = re.sub(
            r"(\bhref\s*=\s*[\"'])([^\"']+)([\"'])",
            update_link,
            original,
            flags=re.IGNORECASE,
        )
        if updated != original:
            page.write_text(updated, encoding="utf-8")
    pdf.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="prepared manual publication tree")
    archive_manual(parser.parse_args().root)
