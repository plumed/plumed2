Dynamic residue-center wrapping regression
==========================================

Known whole physical coordinates are independently shifted by arbitrary lattice
vectors. The oracle uses the known whole coordinates to calculate batch centers
and their expected periodic translations; it does not reconstruct input atoms
or use reference geometry. Rectangular and skew cells, mixed three-/four-atom
residues and an ion, partial/reversed selections and selected-run boundaries
are covered. A second frame changes both geometry and cell.

Metadata fixtures have deliberately unrelated coordinates, and one has a wrong
reference box: output must be identical. Chains, insertion codes, TER boundaries
and nonconsecutive reused identifiers (including a sparse selection) distinguish
residue occurrences. Stable duplicate removal, unchanged unselected atoms,
PAIR count/order/repeated references, virtual COM, overlapping ATOMS/AROUND,
STRIDE and empty selections are checked. A 128-reference case checks the cell
search against the independent direct-distance oracle. Ordinary wrapping stays
per atom; no residue grouping is inferred without the opt-in flag.

Invalid ordinary MOLINFO/atom metadata, virtual selected atoms, GROUPBY conflicts
and both new/legacy PAIR count diagnostics are exercised. No WHOLE or
SMALL_RESIDUES flag, reference cell or whole reference geometry is required.
The committed metadata.pdb also supplies the executable manual example.
