# comparemol

This package is used to compare if two molecules are equivalent.

## Installation

This package requires Python 3.8.

```sh
pip install comparemol
```

## Usage

### Compare if two molecules are equivalent

```py
from comparemol import Mol
# define types and coordinates
types1 = [0, 0]
types2 = [0, 0]
coords1 = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
coords2 = [[0.0, 0.0, 0.0], [1.73205081, 0.0, 0.0]]
# define two molecules
mol1 = Mol(types1, coords1)
mol2 = Mol(types2, coords2)
# compare whether they are equivalent
assert mol1 == mol2
```

### Rotate forces to match coordinates

```py
from comparemol import get_rotation
# define rotation: mol2 -> mol1
r = get_rotation(mol1, mol2)
# rotate forces
force2 = [[1.73205081, 0.0, 0.0]]
force1 = r.apply(force2)
```
