# Solvent Occlusion

Calculate percentage loss of solvent accessible surface area of each residue in a protein upon binding to its receptor, i.e. solvent occlusion. A useful feature for locating "hot spot" residues on PPI interfaces.

## Usage

```python
import solvent_occlusion as so
atoms = so.parse("path/to/pdb/file")
occlusions = so.get_occlusions(atoms)
```

Please also refer to `examples/example.ipynb` for an example.
