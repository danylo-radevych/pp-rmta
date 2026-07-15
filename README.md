# PP-RMTA

Standalone pseudopotential (PP) rigid muffin-tin approximation (RMTA) 
implementation that can be compiled by linking to any installation of 
Quantum ESPRESSO v7.5+, which is used as a DFT engine.

---
Danylo Radevych<sup>1</sup>

1. *Department of Physics and Astronomy,
   George Mason University, Fairfax, VA, USA*



D\. Radevych, T. Shishidou, M. Weinert, E. R. Margine, A. N. Kolmogorov, I. I. Mazin,
*Rigid muffin-tin approximation in plane-wave codes for fast modeling of
phonon-mediated superconductors*,
[npj Comput Mater (2026)](https://doi.org/10.1038/s41524-026-02141-7)

---

## Installation
- Compile [Quantum ESPRESSO (QE)](https://gitlab.com/QEF/q-e.git) in the corresponding `QE_ROOT` folder with:

```
./configure [OPTIONS]
make pw ld1 pp
make install
```

- In the `PPRMTA_ROOT` root folder of the present code, edit the `QE_ROOT` in `make.inc` to link the compiled QE installation, and run

```
make all
```

## Workflow

SCF `$QE_ROOT/bin/pw.x` ---> RMTA `$PPRMTA_ROOT/bin/rmta.x`

## Executable
`rmta.x` in the installation `bin` folder ---> calculate electronic McMillan-Hopfield factors

```
mpirun -n $SLURM_NTASKS $QE_ROOT/bin/pw.x < $prefix.scf.in > $prefix.scf.out
mpirun -n 1 $PPRMTA_ROOT/bin/rmta.x < $prefix.rmta.in > $prefix.rmta.out
```

---

## Pseudopotentials
Hamann's 
[Optimized Norm-Conserving Vanderbilt Pseudopotentials (ONCVPSP)](https://www.pseudo-dojo.org/) 
containing at least `PP_CHI` and `PP_BETA`
blocks are **required**.

**Recommended** ONCVPSP pseudopotentials with explicit `PP_SEMILOCAL` blocks,
regenerated with the [ONCVPSP](https://github.com/oncvpsp/oncvpsp.git) 
version 3.3.1 code based on the input from [PseudoDojo](https://www.pseudo-dojo.org/) 
pseudopotentials, are available in the [oncvpsp-sl](https://github.com/danylo-radevych/oncvpsp-sl.git)
repository. If present, `PP_SEMILOCAL` blocks are used.
Otherwise, semilocal parts are 
recalculated from `chi` functions and `beta` projectors.


---

## RMTA input
| Type         | Variable     | Default   | Description                        |
| :---         | :---         | :---:     | :---                               |
| `CHARACTER`  | `prefix`     | `'pwscf'` | prefix used in the SCF calculation |
| `CHARACTER`  | `outdir`     | `'./'`    | folder where SCF wavefunctions and charge density are stored |
| `LOGICAL`    | `lwrite_dat` | `.false.` | if `.true.`, write spherical potentials and radial functions in `*.dat` files for subsequent plotting with the `plot_lwrite_dat.py` script (not required) |
| `LOGICAL`    | `lhybrid`    | `.false.` | if `.true.`, evaulate energy derivatives of log derivatives as integrals of u^2(r); if `.false.`, evaluate energy derivatives explicitly |
| `LOGICAL`    | `lrmt`       | `.false.` | if `.true.`, read MT radii for each atom, specified in the `rmt(:)` array, from input; if `.false.`, calculate MT radii by dividing nearest-neighbor distances in ratios of the corresponding internal default MT radii, without reading `rmt(:)` |
| `REAL`       | `rmt(:)`     | `-1.0`    | MT radii of each atom in bohr; when `lrmt == .false.`, defaults to MT radii calculated from nearest-neighbor distances and internal default MT radii; if MT radii of any two atoms of the same symmetry type are different, the code stops |

---

<br/><br/>


## Examples
<!-- 
**`$PPRMTA_ROOT/examples/nb-bcc_sl_rmt_default`**: bcc Nb simple metal with 
automatically calculated MT radii, a pseudopotential explicitly
containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`
-->

**`$PPRMTA_ROOT/examples/v3si-a15_sl_rmt_default`**: *A*15 V<sub>3</sub>Si 
compound with automatically calculated MT radii, pseudopotentials
explicitly containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

<!-- 
**`$PPRMTA_ROOT/examples/v-bcc`**: bcc V simple metal with standard pseudopotential from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/v-bcc_sl`**: bcc V simple metal with a pseudopotential explicitly
containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/nb-bcc`**: bcc Nb simple metal with standard pseudopotential from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/nb-bcc_sl`**: bcc Nb simple metal with a pseudopotential explicitly
containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/mo-bcc`**: bcc Mo simple metal with standard pseudopotential from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/mo-bcc_sl`**: bcc Mo simple metal with a pseudopotential explicitly
containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/pd-fcc`**: fcc Pd simple metal with standard pseudopotential from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/pd-fcc_sl`**: fcc Pd simple metal with a pseudopotential explicitly
containing the `PP_SEMILOCAL` block, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/nbn-b1`**: *B*1 NbN compound with standard pseudopotentials from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/nbn-b1_sl`**: *B*1 NbN compound with pseudopotentials explicitly
containing the `PP_SEMILOCAL` blocks, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/v3si-a15`**: *A*15 V<sub>3</sub>Si compound with standard pseudopotentials from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/v3si-a15_sl`**: *A*15 V<sub>3</sub>Si compound with pseudopotentials explicitly
containing the `PP_SEMILOCAL` blocks, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/v3sn-a15`**: *A*15 V<sub>3</sub>Sn compound with standard pseudopotentials from
PseudoDojo, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`

**`$PPRMTA_ROOT/examples/v3sn-a15_sl`**: *A*15 V<sub>3</sub>Sn compound with pseudopotentials explicitly
containing the `PP_SEMILOCAL` blocks, and `lwrite_dat = .true.`; 
plot radial functions with `python plot_lwrite_dat.py`
-->
