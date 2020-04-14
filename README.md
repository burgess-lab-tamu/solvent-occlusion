# Solvent Occlusion

Calculate percentage loss of solvent accessible surface area of each residue in a protein upon binding to its receptor, i.e. solvent occlusion. A useful feature for locating "hot spot" residues on PPI interfaces.

## Installation

Download [the latest release](https://github.com/lyu18/solvent-occlusion/releases/latest), then

```bash
pip install /path/to/the/wheel/file.whl
```

## Usage

```python
import solvent_occlusion as so
atoms = so.parse("path/to/pdb/file")
occlusions = so.get_occlusions(atoms)
```

Please also refer to `examples/example.ipynb` for an example.
