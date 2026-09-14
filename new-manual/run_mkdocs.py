"""Run the MkDocs 2 CLI with its dev5 relative-link export restored."""

import mkdocs
import mkdocs.extensions.relative_urls as relative_urls
from mkdocs.mkdocs import cli, link_to


# MkDocs 2.0.dev5 uses this attribute in its relative_urls extension but does
# not export it from mkdocs.__init__. Remove this compatibility line once the
# pinned release includes the export itself.
if not hasattr(mkdocs, "link_to"):
    mkdocs.link_to = link_to


class _ProtectedMarkdownURL:
    """A non-relative URL marker for Markdown's STX/ETX placeholders."""

    is_relative_url = False


_httpx_url = relative_urls.httpx.URL


def _safe_url(value, *args, **kwargs):
    # Python-Markdown temporarily protects auto-linked email addresses with
    # control characters. MkDocs 2.0.dev5's URL pass sees them too early.
    if "\x02" in value or "\x03" in value:
        return _ProtectedMarkdownURL()
    return _httpx_url(value, *args, **kwargs)


relative_urls.httpx.URL = _safe_url


if __name__ == "__main__":
    cli()
