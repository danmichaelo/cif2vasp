# cif2vasp

This script converts a file from the CIF format to VASP's POSCAR format.
It does not support fractional occupancies. There may also be other 
limitations not yet discovered.

## Usage
  
    $ cif2vasp [-v] [-f] [-e gulp|ase|cctbx] filename.cif

## Backends

The first version of this script was inspired by the procedure found at
http://encina.northwestern.edu/index.php/Cif_to_VASP, but 
I added a better CIF parsing routine using the PyCifRW library.
While this version of the script (to be referred to as the GULP backend) 
still works great, it has very many dependencies.
In order to simplify things, I first looked into CCTBX, and later ASE,
both powerful Python packages, with ASE being the simplest one to install
and use.

### ASE

The excellent ASE package, available from https://wiki.fysik.dtu.dk/ase/,
almost eliminates the need for this script, since converting a file actually is done with just three lines of code:

    from ase import io
    atoms = io.read(jobname+'.cif')
    atoms.write('POSCAR', format = 'vasp')

Using other backends is recommended if you run into problems with ASE's 
CIF parsing, or want direct coordinates directly (ASE writes cartesian).

### GULP 

This backend will give slightly more feedback during the conversion, and
will provide a warning if fractional occupancies are found, but the 
numerical noise is about 1e-10, compared to 1e-16 for ASE because of the
intermediate XYZ file with only 10 digits.
The GULP process will write direct coordinates.

Requirements:

 *  GULP: 
    https://www.ivec.org/gulp/
    (tested with version 3.4.9, other versions may work as well)

 *  AFLOW (aconvasp): 
    http://nietzsche.mems.duke.edu/aflow.html
    (tested with version 30060, other versions may work as well)

 *  PyCifRW 3.3 
    http://pycifrw.berlios.de/ 
    (older versions are not likely to work):

Process:

 1.  `prepareGulpInput`: Read CIF file using the PyCifRW library and write GULP in file. 
 3.  `runGulp`: Run GULP to generate all unit cell coordinates based on the spacegroup symmetry operations.
 4.  `runConvasp`: Write XYZ file from GULP out file, then use ACONVASP to convert the cartesian coordinates 
      in the XYZ file into fractional coordinates.

I also experimented with ATAT (ezvasp) from http://www.its.caltech.edu/~avdw/atat/ .
You may try modifying the script to run `cif2vaspUsingGULP(jobname, ezvasp = True)`

### CCTBX

I looked into using the Computational Crystallography Toolbox (CCTBX), 
which is a powerful Python package available from http://cctbx.sourceforge.net.
I did however run into some problems getting CCTBX to play nicely with my current Python 
installation. The CCTBX support is therefore experimental.

