"""Post-process a MkDocs 2 build for compatibility with existing manual URLs."""

from pathlib import Path


def redirect_html(target):
    return f"""<!doctype html>
<html lang=\"en\">
  <head>
    <meta charset=\"utf-8\">
    <meta http-equiv=\"refresh\" content=\"0; url={target}\">
    <link rel=\"canonical\" href=\"{target}\">
    <title>Redirecting…</title>
  </head>
  <body><p><a href=\"{target}\">Continue to the requested page</a></p></body>
</html>
"""


def main():
    docs = Path("docs")
    site = Path("site")
    for source in docs.rglob("*.md"):
        stem = source.stem
        if source.name.lower() in ("index.md", "readme.md") or stem == stem.lower():
            continue
        relative_parent = source.relative_to(docs).parent
        old_dir = site / relative_parent / stem
        new_dir = site / relative_parent / stem.lower()
        if old_dir == new_dir or not new_dir.exists():
            continue
        # On case-insensitive filesystems both spellings are the same directory.
        if old_dir.exists() and old_dir.samefile(new_dir):
            continue
        old_dir.mkdir(parents=True, exist_ok=True)
        depth = len(old_dir.relative_to(site).parts)
        target = "../" * depth + "/".join(new_dir.relative_to(site).parts) + "/"
        (old_dir / "index.html").write_text(redirect_html(target))


if __name__ == "__main__":
    main()
