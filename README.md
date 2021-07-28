# Quickstart

## Installation

Install this module with pip

```sh
pip install -U git+https://github.com/mtitze/touschek.git
```

## Usage

```python
from touschek import optics_tools

lattice_filename = 'name_of_your_lattice'

# Set some beam input parameters. Note that EX and EY may be changed later, when the natural parameters
# are computed:
beam = {'PARTICLE': 'ELECTRON', 'ENERGY': 1.7, 'EX': 7.54543712e-09, 'EY': 1.2e-10, 'NPART': 5000000000.0, 'SIGE': 0.00069541, 'SIGT': 0.0029196, 'radiate': True}

# Create optics class with given lattice and beam parameters:
opt = optics_tools.optics(lattice_filename, beam_params=beam, show_init=False, verbose=True)

# Compute Touschek-lifetime, here assuming a certain machine symmetry for speed increase and using
# fast routine (precise == False). The optional coupling parameter coupling_y can be used to enlarge EY by coupling_y*EX:
touschek_results = opt.touschek_lifetime(precise=False, symmetry=8, coupling_y=0.02)
```

## Further reading

https://touschek.readthedocs.io/en/latest/index.html
