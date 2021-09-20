# Examples

Here you can find a set of examples showing how the code can be used. They are separated
in the directories `mpi` (multi-core examples) and `nonmpi` (single-core examples).
Inside each directory you will find examples of synthesis and inversions.

The script `syn_obs.py` will generate all observations needed for the following scripts.

## nonmpi/syn

#### `example_prog.py`
Programmatic definition of a chromosphere and how to get the synthetic Stokes profiles when changing some of the parameters.

#### `example_single_syn.py`
Basic example showing how to do a single-pixel synthesis.

#### `example_nonmpi_syn1d.py`
Basic example showing how to do a single-pixel synthesis but using iterators. This can be used to do
multi-pixel synthesis in a serial machine.

#### `example_coordinates.py`
Basic example showing how different coordinates systems can be used.

#### `caii_syn.py`
Basic example showing how to synthesize Ca II 8542 A in NLTE.

#### `check_boundary.py`
Basic example showing the effect of changing the boundary condition in chromospheres.

## nonmpi/inv

#### `example_single_inv.py`
Basic example showing how to do the inversion of a single pixel.

#### `example_nonmpi_inv1d.py`
Basic example showing how to do the inversion of a single pixel (although it can deal with many pixels) using the machinery
for the inversion of many pixels.

#### `example_algorithm.py`
Basic example showing how to use a different algorithm to do the inversion.

#### `example_coupled.py`
Basic example showing how to carry out inversions in He I 10830 and D3 lines simultaneously.

## mpi/syn

#### `example_mpi_synh5.py`
Basic example showing how to carry out synthesis in parallel.

## mpi/inv

#### `example_mpi_invh5.py`
Basic example showing how to carry out inversions in parallel.

#### `example_mpi_invh5_mask.py`
Basic example showing how to carry out inversions in parallel and masking some pixels.