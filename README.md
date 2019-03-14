omm-top: Gromacs topology parser for OpenMM
===========================================

```python
import ommtop
from simtk.openmm import app
from simtk.unit import nanometer

gro = GromacsGroFile('input.gro')
top = ommtop.GromacsTopFile(
    'input.top',
    periodicBoxVectors=gro.getPeriodicBoxVectors(),
    includeDir='/usr/local/gromacs/share/gromacs/top',
    defines={'FLEXIBLE': ''},
)
system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*nanometer,
        constraints=app.HBonds)
```

The parser aims at being mostly compatible with the one provided in OpenMM,
while implementing more Gromacs potentials and features.
