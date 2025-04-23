The variables within this module can be used to detect the presence of secondary structure elements such as alpha helices or parallel and antiparallel beta sheets in a protein.

Notice that the [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md) are implemented as shortcuts.  In other words, the only non shortcut
action in this module is the [SECONDARY_STRUCTURE_RMSD](SECONDARY_STRUCTURE_RMSD.md) that calculates the RMSD distance between each segement of residue and the idealised version of 
secondary structure element.  This structure makes it easy to add additional CVs that work similarly to the [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md) 
either directly in the input file or by adding a new shortcut action.
