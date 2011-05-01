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

 *  AFLOW (aconvasp): 
    http://nietzsche.mems.duke.edu/aflow.html

 *  ATAT (ezvasp): 
    http://www.its.caltech.edu/~avdw/atat/

 *  PyCifRW 3.3: 
    http://pycifrw.berlios.de/ 

I also tried to make a version of cif2vasp that made use of the
Computational Crystallography Toolbox (CCTBX), but I ran into some
problems with using CCTBX with my Python installation. 
CCTBX is available from
  http://cctbx.sourceforge.net
Have a look at the `cif2vaspUsingCCTBX` function if you want to play with it.

## Usage
  
    $ cif2vasp filename.cif

