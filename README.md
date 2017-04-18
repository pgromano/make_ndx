# make_ndx

Python API which uses the parmed toolkit to parse index values for molecular topology files.

### Installation
Clone the package and change directory into the `make_ndx` directory.

    python -m pip install -e . --user

### Usage
While the `Parser` object has additional attributes to parse through atom, residue, and chain sequences; the primary purpose is to find the index `where` atomic sites are located.

    In [1]: from make_ndx import Parser

    In [2]: p = Parser('molecule.pdb')
    In [3]: p.where('C5')
    Out[3]: array([ 8, 25])

    In [4]: p.where('C5', resid=1)
    Out[4]: array([ 25])
