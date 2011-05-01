# cif2vasp

This script converts a file from the CIF format to VASP's POSCAR format.
It does not work with fractional coordinates. There may also be other 
limitations not yet discovered.

This script is inspired by the procedure found at
http://encina.northwestern.edu/index.php/Cif_to_VASP

## Installation

This script is NOT self-contained. It makes use of external 
programs and python libraries. Please make sure you have installed
recent versions of the following software:

 *  GULP: 
    https://www.ivec.org/gulp/
    (tested with version 3.4.9, other versions may work as well)

 *  AFLOW (aconvasp): 
    http://nietzsche.mems.duke.edu/aflow.html
    (tested with version 30060, other versions may work as well)

 *  PyCifRW 3.3 
    http://pycifrw.berlios.de/ 
    (older versions are not likely to work):

## Optional packages 

It would be preferable to get rid of the dependency on GULP and AFLOW, 
especially GULP which is a quite big package that takes some time to compile.

One alternative I've looked into is the Computational Crystallography Toolbox (CCTBX), 
which is a very powerful Python package available from http://cctbx.sourceforge.net.
I ran into some problems implementing CCTBX with my current Python installation, but
an experimental method `cif2vaspUsingCCTBX` is available as an alternative to
`cif2vaspUsingGULP`.

I also experimented with ATAT (ezvasp) from http://www.its.caltech.edu/~avdw/atat/ .
You may try modifying the script to run `cif2vaspUsingGULP(jobname, ezvasp = True)`

## Usage
  
    $ cif2vasp filename.cif

## Code outline

 1.  `prepareGulpInput`: Read CIF file using the PyCifRW library and write GULP in file. 
 3.  `runGulp`: Run GULP to generate all unit cell coordinates based on the spacegroup symmetry operations.
 4.  `runConvasp`: Write XYZ file from GULP out file, then use ACONVASP to convert the cartesian coordinates 
      in the XYZ file into fractional coordinates.
