# Installing

Install the documentation dependencies in a virtual environment:

```sh
python -m pip install -r requirements.txt
```

The manual uses the MkDocs 2 pre-release pinned in `requirements.txt`. Its
templates and browser-side functionality are kept in this repository rather
than supplied by an external theme.

# Building 

You then build the website locally by doing

```sh
make
```

To see the website you do

```sh
python run_mkdocs.py serve
```

Run `make` before starting the preview server, because the Markdown pages and search index
are generated from the compiled PLUMED executable.
