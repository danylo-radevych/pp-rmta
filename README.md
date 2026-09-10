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

- In the `PPRMTA_ROOT` root folder of the present code, edit the `QE_ROOT` in `make.inc` to link the compiled QE installation, along with optional additional compiler flags in `OPTFLAGS`, and run

```
make all
```

## Workflow

SCF `$QE_ROOT/bin/pw.x` ---> [NSCF `$QE_ROOT/bin/pw.x`] ---> RMTA `$PPRMTA_ROOT/bin/rmta.x`

## Executable
`rmta.x` in the installation `bin` folder ---> calculate electronic McMillan-Hopfield factors

```
 mpirun -n $SLURM_NTASKS $QE_ROOT/bin/pw.x < $prefix.scf.in > $prefix.scf.out
[mpirun -n $SLURM_NTASKS $QE_ROOT/bin/pw.x < $prefix.nscf.in > $prefix.nscf.out]
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
| `CHARACTER`  | `formulation`| `'default'` | if `'nodeless'`, calculate partial DOS based on continuous wavefunctions; if `'derivative'`, calculate partial DOS based on continuous derivatives of the wavefunctions; if `'default'`, automatically activate `'derivative'` when `'nodeless'` is close to radial nodes (recommended) |
| `CHARACTER`  | `prefix`     | `'pwscf'` | prefix used in the SCF calculation |
| `CHARACTER`  | `outdir`     | `'./'`    | folder where SCF wavefunctions and charge density are stored |
| `CHARACTER`  | `rmt_method` | `'touching'`| see `lrmt`: if `'default'`, set default MT radii from the table; if `'pseudo'`, set default MT radii as pseudopotential cutoff radii; if `'neighbor'`, divide nearest-neighbor distances into ratios of the default MT radii; if `'touching'`, starting from `'neighbor'`, enforce touching spheres (recommended); if `'pseudoneighbor'`, divide nearest-neighbor distances into ratios of the pseudopotential cutoff radii; if `'pseudotouching'`, starting from `'pseudoneighbor'`, enforce touching spheres|
| `LOGICAL`    | `lrmt`       | `.false.` | if `.true.`, read MT radii for each atom, specified in the `rmt(:)` array, from input; if `.false.`, calculate MT radii automatically with one of the methods specified in `rmt_method`, without reading `rmt(:)` |
| `LOGICAL`    | `ltetra`     | `.true.`  | if `.true.`, use tetrhedron method for partial DOS integration; if `.false.`, use delta-function smearing set by `ngauss` and `degauss` |
| `LOGICAL`    | `lwrite_dat` | `.false.` | if `.true.`, write spherical potentials and radial functions in `*.dat` files for subsequent plotting with the `plot_lwrite_dat.py` script (not required) |
| `LOGICAL`    | `lhybrid`    | `.true.`  | if `.true.`, evaulate radial integrals of u^2(r, e) explicitly; if `.false.`, evaluate them through the derivatives of u(r, e) |
| `INTEGER`    | `ngauss`     | `-99`     | delta-function smearing for partial DOS integration: `-99` Fermi-Dirac; `0` Gauss; `1` MP; see `ltetra`, `degauss` |
| `INTEGER`    | `igp_scale`  | `3`       | integer scaler for the number of Gauss points used in spherical integration |
| `REAL`       | `degauss`    | `0.001`   | smearing degauss value in Ry; see `ltetra`, `ngauss` |
| `REAL`       | `rmt(:)`     | `-1.0`    | MT radii of each atom in bohr; when `lrmt == .false.`, defaults to MT radii calculated with a method in `rmt_method`; if MT radii of any two atoms of the same symmetry type are different or atomic spheres overlap, the code stops |


---

<br/><br/>


## Examples


**`examples/nb-bcc_sl_rmt_touching`**: bcc Nb simple metal with automatically calculated MT radii and a pseudopotential explicitly containing the `PP_SEMILOCAL` block

**`examples/v3si-a15_sl_rmt_touching`**: *A*15 V<sub>3</sub>Si compound with automatically calculated MT radii and pseudopotentials explicitly containing the `PP_SEMILOCAL` blocks

**`examples/v3si-a15_nscf_sl_rmt_touching`**: *A*15 V<sub>3</sub>Si compound with automatically calculated MT radii and pseudopotentials explicitly containing the `PP_SEMILOCAL` blocks. SCF is performed on a coarse *k*-point grid, and optional NSCF step is used to get a finer *k*-point grid for the subsequent RMTA calculation.





