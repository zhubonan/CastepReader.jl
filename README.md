# CastepReader.jl

Julia package for reading CASTEP input/output files.

Implemented readers:

- `ome_bin`
- `cst_ome`
- `bands`
- `cell`
- `param`
- `pdos_bin`

Working but still being improved:

- `castep_bin`/`check`


## Notes

- Units are not handled for binary data from check/castep_bin file - these are the internal atomic units 
  - Energies are likely to be in Hartree
  - Length units are likey to be in Bohr
  - Lattice vectors are stored as matrices of row vectors rather than the column vectors

- The unit system may need to be reworked - returning data with units may not always be desirable (for sake type dispatch)

- More quantities from the `castep_bin/check` file can be read if the scehma is known.
- Units of the wavefunction/charge density needs to be determined
