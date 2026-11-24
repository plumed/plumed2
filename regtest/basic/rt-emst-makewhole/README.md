Regression for issue #1496: MOLINFO WHOLE must build the reconstruction tree
and pair each parent with the following entry in the tree-index list.

Four atoms form a line at x = 1, 7, 3, 5 Angstrom. Their minimum spanning tree
visits atoms 1, 3, 4, 2, with edges of length 2 Angstrom. Sequential atom-list
reconstruction is incorrect in the 10 Angstrom box.

The first two frames split this line across periodic images; the second frame
also translates every input atom by one box length, checking tree reuse and
anchor preservation. The third frame is already whole. No WHOLEMOLECULES action
runs before the tested GYRATION, CENTER and RMSD actions. A reordered subset
checks local indexes independently of absolute atom indexes.

The whole line has center x=4 and radius of gyration sqrt(5), with zero RMSD.
The translated frame has center x=14. The subset has radius sqrt(8/3); its center
is -7 in the first frame and 3 in the other two, because its first atom anchors
the image. NOPBC controls retain the split coordinates.
