# cif2vasp

This script creates INCAR POSCAR POTCAR and KPOINTS from the 
filename.cif and potentials located in [...]

  $ cif2vasp filename.cif

## Installation

This script is NOT self-contained, since it makes use of external 
programs and python libraries. Please make sure you have installed
recent versions of the following software:

  GULP: 
  https://www.ivec.org/gulp/

  AFLOW (aconvasp): 
  http://nietzsche.mems.duke.edu/aflow.html

  ATAT (ezvasp): 
  http://www.its.caltech.edu/~avdw/atat/

  PyCifRW 3.3: 
  http://pycifrw.berlios.de/ 

Optional:

  Computational Crystallography Toolbox (CCTBX)
  http://cctbx.sourceforge.net

## Usage

Create your INCAR file first. Make sure that the file has '[INCAR]' 
at the beginning, KPPRA, DOGGA, KSCHEME, and other relevant ATAT 
input tags at the end.
Also, make sure it has a [POSCAR] appendage as the last line, with 
no newlines afterward.

E.g.
[INCAR]
-Incar stuff that you want

KSCHEME = MP
DOGGA
KPPRA = 10000
[POSCAR]

