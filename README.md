# mof_vac_level

This is an installable version of the code proposed by Fumanal et al. [1], building on top of work from Walsh and co-workers [2], as implemented in the [MacroDensity Python package](https://github.com/WMD-group/MacroDensity).

⚠️ **Not well tested**   ⚠️


## Changes

- object-oriented API
- some parts of the code are vectorized, just-in-time compiled or multiprocessed
- progress bar

## Installation

```
pip install git+git@github.com:kjappelbaum/mof_vacuum_level.git
```

## Usage

```(python)
from mof_vac_level import MOFVacLevel
mvl = MOFVacLevel('aiida-ELECTRON_DENSITY-1_0.cube')
mvl.get_vacuum_potential(res=0.4)
```

## References
[1] [Fumanal, M.; Capano, G.; Barthel, S.; Smit, B.; Tavernelli, I. Energy-Based Descriptors for Photo-Catalytically Active Metal–Organic Framework Discovery. J. Mater. Chem. A 2020, 8 (8), 4473–4482.](https://doi.org/10.1039/C9TA13506E)

[2] [Butler, K. T.; Hendon, C. H.; Walsh, A. Electronic Chemical Potentials of Porous Metal–Organic Frameworks. J. Am. Chem. Soc. 2014, 136 (7), 2703–2706.](https://doi.org/10.1021/ja4110073)
