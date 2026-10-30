This directory is **not** used to build the user manual.

The user manual is built with MkDocs from `new-manual/` (`make doc`), and is
published and installed as `user-doc/html`, which is why that path still
appears in links throughout the code.

What is left here comes from the old Doxygen manual and is not rendered anywhere:

- `tutorials/` and `figs/`: old tutorials (Doxygen `.txt` pages) with their inputs
  and figures. Current tutorials live at https://www.plumed-tutorials.org
- `Performances.md`: hints on improving performance, not yet migrated to the new manual.

New documentation should not be added here.
