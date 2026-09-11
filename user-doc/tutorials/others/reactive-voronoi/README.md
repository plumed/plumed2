# Reactive Soft-Voronoi runnable examples

This directory provides two small, full-pair `plumed driver` checks for the
three Reactive Soft-Voronoi Actions. The fixtures are designed to make atom
numbering, `CENTERS`, `ASSIGNED`, and `REFERENCE` auditable. They are teaching
and implementation checks, not condensed-phase production models.

## Files

| File | Purpose |
| --- | --- |
| `two-water.xyz` | Two neutral water molecules in O-H-H order. |
| `water.plumed` | Broadcast `REFERENCE=2` and evaluate water defect activity and separation. |
| `glycine-two-water.xyz` | Zwitterionic glycine with one complete water near N and one near a carboxyl O. |
| `glycine.plumed` | Evaluate the paper-consistent glycine protonation and defect-separation coordinates. |

## Prerequisite

Use a PLUMED build that contains `VORONOI_COORDINATION`, `VORONOI_DISTANCE`,
and `VORONOI_POSITION`. These Actions belong to the default `colvar` module,
so the supplied inputs use them directly and contain no `LOAD` line. Ensure
that `plumed driver` resolves to the installation being tested before running
the examples.

When reviewing the source before it is available in an installed PLUMED
version, an out-of-tree runtime plugin can still be used as a development
fallback. Build it with the same PLUMED executable and compiler/ABI used by the
target program and add the corresponding `LOAD` line only to the local test
input. Do not commit copied source files, shared libraries, or local `LOAD`
paths to this directory.

## Example 1: two neutral waters

The XYZ atom order is:

| Atoms | Selection | Internal REFERENCE |
| --- | --- | --- |
| `1,4` | `WaterO`, `CENTERS` | `(2,2)` |
| `2,3,5,6` | `WaterH`, `ASSIGNED` | - |

The scalar input `REFERENCE=2` is broadcast to both water O centers. It is
shorthand for `(2,2)`, not a one-center calculation.

Run:

```bash
plumed driver --plumed water.plumed --ixyz two-water.xyz \
  --box 12.028,12.028,12.028 --length-units A
```

The single neutral frame gives:

```text
#! FIELDS time sa st
 0.000000 0.000000 0.000000
```

Soft values become nonzero on proton-transfer and H₃O⁺/OH⁻ frames. This neutral
fixture checks indexing and the zero-defect limit only.

`WATER_DERIVATIVES` contains 28 lines for the two Actions on this six-atom
frame.

## Example 2: zwitterionic glycine and two waters

The translated, renumbered fixture was extracted from an optimized 128-water
zwitterionic geometry. It retains the closest complete water to glycine N and
the closest complete water to a carboxyl O.

| Atoms | Selection and role | REFERENCE |
| --- | --- | --- |
| `1,4` | `WaterO`, `CENTERS` | `2,2` |
| `2,3,5,6` | `WaterH`, `ASSIGNED` | - |
| `7` | `GlyN`, `CENTER` | `2` |
| `8,9,10` | Three N-H atoms, `ASSIGNED` | - |
| `11-14` | Glycine C/C-H framework, not selected | - |
| `15,16` | `GlyO1,GlyO2`, `CENTERS` | `0.5,0.5` |

`CENTERS` expands to `(1,4,7,15,16)`. The reference vector
`(2,2,2,0.5,0.5)` defines the neutral [N] chemical origin even though the
coordinate fixture is zwitterionic [Z]. In an ideal [Z] frame, the five
defects are approximately `(0,0,+1,-0.5,-0.5)`. The reference sum is seven,
equal to the four water H and three N-H atoms in `AllH`.

Run:

```bash
plumed driver --plumed glycine.plumed --ixyz glycine-two-water.xyz \
  --box 15.813,15.813,15.813 --length-units A
```

The supplied frame gives:

```text
#! FIELDS time sp sd_water sd_internal sd
 0.000000 0.007332 0.035428 2.978850 3.014278
```

Thus `sp` is near zero and `sd` is near 3 Å, as expected for [Z]. The small
nonzero deviations reflect the finite KAPPA soft assignment rather than an
indexing error. This two-water excerpt cannot reproduce the condensed-phase
free-energy surface reported for solvated glycine.

`GLYCINE_DERIVATIVES` contains 46 lines for the three differentiable
Reactive Soft-Voronoi Actions on this 16-atom frame.

## Choosing KAPPA and NLIST

- A larger KAPPA approaches a hard nearest-center assignment and can sharpen
  printed state labels. A smaller KAPPA generally gives smoother CV and force
  changes for biased sampling, but can blur distinct basins if it is too small.
- Both examples use `KAPPA=5`. Values such as 5 or 100 are system-specific
  choices, not universal defaults; scan representative reactant, transition,
  product, and host-switching frames before biasing.
- These fixtures intentionally omit `NLIST`, so every CENTER-ASSIGNED pair is
  included and there is no hidden cutoff. If `NLIST` is enabled, provide both
  `NL_CUTOFF` and `NL_STRIDE`, then converge values and derivatives against the
  exact full-pair calculation.

The glycine definitions follow [Intramolecular and Water Mediated Tautomerism
of Solvated Glycine](https://doi.org/10.1021/acs.jcim.4c00273). An additional
application is described in [Modulation of Electric Field and Interface on
Competitive Reaction Mechanisms](https://doi.org/10.1021/acs.jctc.5c00705).
