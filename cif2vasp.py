#! /usr/bin/env python
#
# Inspired by: http://encina.northwestern.edu/index.php/Cif_to_VASP
#
# See README.md
#

import sys, os, re
import subprocess
from optparse import OptionParser

#from cctbx import uctbx, sgtbx, crystal
#from cctbx import xray
#from cctbx import crystal
#from cctbx.array_family import flex

def usage():
    print "Usage: cif2vasp.py [-v] [-f] [-e gulp|ase|cctbx] filename.cif "
    print "Options"
    print " -v, --verbose  increase verbosity"
    print " -f             try to identify fractional numbers (special positions) "
    print "                from coordinates with less than six decimals."
    print " -e, --engine   the backend to use for the conversion (gulp, ase or cctbx)"

def readCifFile(cifFile):
    from CifFile import CifFile

    if not os.path.exists(cifFile):
        raise IOError("CIF file '%s' was not found!" % (cifFile))
    
    cf = CifFile(cifFile)
    print "------------------------------------------------------------------"
    if len(cf) != 1:
        raise StandardError("The cif file contains %i data blocks, while one was expected")
        # A cif file can contain several "datablocks" that each start
        # with "data_".
    
    cb = cf[cf.keys()[0]]                               # open the first block
    AA = float(re.match('([0-9.]*)',cb['_cell_length_a']).group(0))
    BB = float(re.match('([0-9.]*)',cb['_cell_length_b']).group(0))
    CC = float(re.match('([0-9.]*)',cb['_cell_length_c']).group(0))
    alpha = float(cb['_cell_angle_alpha'])
    beta = float(cb['_cell_angle_beta'])
    gamma = float(cb['_cell_angle_gamma'])
    SG = int(cb['_symmetry_Int_Tables_number'])              # spacegroup
  
    atomTypes = []
    atoms = ''
    fracOccFound = False
    firstAtom = True
    atoms = []
    for atom in cb.GetLoop('_atom_site_label'):
        atomKeys = dir(atom)
        if '_atom_site_type_symbol' in atomKeys:
            m = re.match('[a-z]*',atom._atom_site_type_symbol,re.I)
            atomType = m.group(0)
        else:
            m = re.match('[a-z]*',atom._atom_site_label,re.I)
            atomType = m.group(0)
        
        atomLabel = atom._atom_site_label

        if '_atom_site_occupancy' in atomKeys:
            occ = float(atom._atom_site_occupancy)
            if not occ == 1.0:
                if not fracOccFound: 
                    print " "
                print "  WARNING: Fractional occupancy (" + str(occ) +") " \
                    + "found for atom of type " + atomType + "."                
                fracOccFound = True
        else:
            occ = 1.0
        
        # Some crystal structures obtained by neutron diffraction use D for H:
        if atomType == 'D':
            atomType = 'H'
            atomLabel.replace('H','D')
        
        if '_atom_site_symmetry_multiplicity' in atomKeys and '_atom_site_Wyckoff_symbol' in atomKeys:
            atomTypes.append(atomType+' at '+atom._atom_site_symmetry_multiplicity+atom._atom_site_Wyckoff_symbol)
        else:
            atomTypes.append(atomType)
        
        atomPos = [atom._atom_site_fract_x, atom._atom_site_fract_y, atom._atom_site_fract_z]
        for p in atomPos:
            pp = p.split(".")
            if len(pp) is 2:
                decimals = p.split(".")[1]
                if len(decimals) > 3 and len(decimals) < 6 and decimals[0] == decimals[1] and decimals[-1] != "0":
                    print "\n  ---------------------\n"\
                        + "  Warning: If the fractional coordinate "+p+" is a recurring decimal, such as 1/3,\n" \
                        + "    then it is necessary to specify this value to six decimal places to be sure of \n" \
                        + "    it being recognised correctly as a spcecial position.\n  ------------------" 
		
		
		# The coordinates of the atom (_atom_site_fract_x/y/z) may have 
		# a last digit in parenthesis, like "0.6636(7)". Therefore we
		# extract the part consisting of only digits and a decimal separator:
		p = re.compile('[0-9.]*');
        atomX = float(p.match(atom._atom_site_fract_x).group())
        atomY = float(p.match(atom._atom_site_fract_y).group())
        atomZ = float(p.match(atom._atom_site_fract_z).group())
        
        #atoms += "%s %f %f %f %f %f\n" % (atomType, atomX, atomY, atomZ, 0.0, occ)
        atoms.append({'label': atomLabel, 'type': atomType, 'pos': (atomX,atomY,atomZ) })
        firstAtom = False

    if fracOccFound: 
        print " "
        print "ERROR: Fractional occupancies are not currently supported.\n"
        exit()
    
    print "  Atom types: " + ', '.join(atomTypes)
    
    return {'spacegroup': SG, 'unit_cell': [AA,BB,CC,alpha,beta,gamma], 'scatterers': atoms}
    
def cif2vaspUsingCCTBX(jobname, ezvasp = False):
    
    #print "WARNING: This script does NOT work with files that have fractional occupancies."

    # os.mkdir(jobname)
    
    cif = readCifFile(jobname+'.cif')

    print "CIF file read successfully:"

    unit_cell = uctbx.unit_cell(cif['unit_cell'])
    space_group_info = sgtbx.space_group_info(symbol=cif['spacegroup'])
    crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,space_group_info=space_group_info)
    crystal_symmetry.show_summary()
    
    print " "
    #print "  Space group:",SG
    #print "  a=%s, b=%s, c=%s, alpha=%s, beta=%s, gamma=%s" % (AA,BB,CC,alpha,beta,gamma)
    
    #print cif['scatterers']
    
    scatterers = flex.xray_scatterer()
    for s in cif['scatterers']:
        scatterers.append(xray.scatterer(label=s['label'], site=s['pos']))

    print
    print "--------------- icsd_structure ---------------"
    print
    icsd_structure = xray.structure(crystal_symmetry=crystal_symmetry, scatterers=scatterers)
    icsd_structure.show_summary().show_scatterers()
    #print 
    #icsd_structure.show_distances(distance_cutoff=2.5)
  
    print
    print "--------------- primitive_structure ---------------"
    print
    primitive_structure = icsd_structure.primitive_setting()
    primitive_structure.show_summary().show_scatterers()

    print
    print "--------------- p1_structure ---------------"
    print
    p1_structure = primitive_structure.expand_to_p1()
    p1_structure.show_summary().show_scatterers()
    print
    print "OK"    
   
# Requires existing gulp out file for unit cell
def xyz2vaspUsingGULP(jobname, ezvasp = False, verbose = False, auto_fractions = False):
    
    # convasp will convert your atom positions into fractional 
    # and produce a file that looks just like a VASP POSCAR:    
    runConvasp(
        jobName = jobname,
        gulpOutputFile = jobname + '.gulp.out',
        xyzFile = jobname+'.xyz',
        ezvaspStyle = ezvasp
    )
    
    if ezvasp:
        prepareEzvaspInput(jobname)
        sys.stdout.write("Creating INCAR, POSCAR, POTCAR, KPOINTS using ezvasp... ")    
        sys.stdout.flush()
        os.system("ezvasp -n vasp.in")
        sys.stdout.write("done\n")
    else:
        os.rename(jobname+'.convasp.out','POSCAR')
    
    
def cif2vaspUsingGULP(jobname, ezvasp = False, verbose = False, auto_fractions = False):

    prepareGulpInput(
        cifFile = jobname + '.cif', 
        gulpFile = jobname + '.gulp.in',
        jobName = jobname,
        verbose = verbose,
        auto_fractions = auto_fractions
    )
    runGulp(jobname, verbose)
    
    # convasp will convert your atom positions into fractional 
    # and produce a file that looks just like a VASP POSCAR:    
    runConvasp(
        jobName = jobname,
        gulpOutputFile = jobname + '.gulp.out',
        xyzFile = jobname+'.xyz',
        ezvaspStyle = ezvasp,
        verbose = verbose
    )
    
    if ezvasp:
        prepareEzvaspInput(jobname)
        sys.stdout.write("Creating INCAR, POSCAR, POTCAR, KPOINTS using ezvasp... ")    
        sys.stdout.flush()
        os.system("ezvasp -n vasp.in")
        sys.stdout.write("done\n")
    else:
        os.rename(jobname+'.convasp.out','POSCAR')

def cif2vaspUsingASE(jobname):
    """
    ASE seems to do an excellent job with reading cif's.
    It will write out the coordinates in cartesian coordinates.
    """
    from ase import io
    atoms = io.read(jobname+'.cif')
    atoms.write('POSCAR', format = 'vasp')


def prepareEzvaspInput(jobname):

    ezvaspIn = open('vasp.in','w')

    incar = file('INCAR')
    while True:
        line = incar.readline()
        if len(line) == 0:
            break
        ezvaspIn.write(line) # notice comma
    incar.close()

    ezvaspIn.write('\n[POSCAR]\n') # notice comma

    poscar = file(jobname + '.convasp.out')
    lineNo = 0
    while True:
        lineNo += 1
        line = poscar.readline()
        if len(line) == 0:
            break
        if lineNo != 6: # ezvasp is not interested in atom counts line
            ezvaspIn.write(line)
    poscar.close()

    ezvaspIn.close()


def findLinesContaining(lines, str):
    enumLines = enumerate(lines)
    return [k for k, v in enumLines if str in v]

# === CONVASP ======================================================================

def runConvasp(gulpOutputFile, xyzFile, jobName, ezvaspStyle = True, verbose = False): 
    """
    Example file:
        mgh2
    1 <-- scaling factor (leave as one)
            4.516800    0.000000    0.000000
            0.000000    4.516800    0.000000
            0.000000    0.000000    3.020500 <-- Cartesian lattice vectors (Angstrom)
    2 4 <-- number of atoms of each type (2 Mg atoms, 4 H atoms)
    Cartesian <-- coordinate style
            0.000000000         0.000000000         0.000000000 Mg
            2.258400000         2.258400000         1.510250000 Mg
             1.382140800         1.382140800         0.000000000 H
             3.134659200         3.134659200         0.000000000 H
             3.640540800         0.876259200         1.510250000 H
             0.876259200         3.640540800         1.510250000 H
    """
    
    xyzfile = open(xyzFile,'r')         #read the GULP output files
    outfile = open(gulpOutputFile,'r')
    XYZ = xyzfile.readlines()
    OUT = outfile.readlines()
    
    # Find cartesian lattice vectors in Gulp Output file:
    lineNo = findLinesContaining(OUT, 'Cartesian lattice vectors (Angstroms)')[0];
    latVec = [
        OUT[lineNo + 2][:-1].split(),
        OUT[lineNo + 3][:-1].split(),
        OUT[lineNo + 4][:-1].split()
    ] # :-1 removes \n
    
    # Find cartesian coordinates in xyz file:
    atoms = []
    atomTypes = []
    atomCounts = []
    atomCount = 0
    for lineNo in range(2, len(XYZ)):
        line = XYZ[lineNo].split()
        atomType = line[0]
        
        if len(line) == 4:  # A valid coordinate line should have length 4
            if ezvaspStyle:
                # Include the element name if we are to use EzVasp ...
                atoms.append("%s %s %s %s" % (line[1],line[2],line[3],atomType))
            else:
                # ... or drop it if we are to use Vasp directly:
                atoms.append("%s %s %s" % (line[1],line[2],line[3]))
        
        if len(atomTypes) == 0:
            atomTypes.append(atomType)
        else:
            if not atomTypes[len(atomTypes)-1] == atomType: # new atom type
                atomTypes.append(atomType)
                atomCounts.append(str(atomCount))
                atomCount = 0
        atomCount += 1
    atomCounts.append(str(atomCount)) # store the count of the last atom type
    
    if verbose:
        for i in range(len(atomTypes)):
            print "  found %s atoms of type %s" % (atomCounts[i],atomTypes[i])
    
    convaspInput = [jobname,'1']
    convaspInput.extend([' '.join(v) for v in latVec])
    convaspInput.extend([' '.join(atomCounts),'Cartesian'])
    convaspInput.extend(atoms)
    stdin = '\n'.join(convaspInput)

    if verbose:
        print "Converting to POSCAR format using aconvasp... "
    p = subprocess.Popen(['aconvasp','--direct'], 
            stdin = subprocess.PIPE, 
            stdout = file(jobname+'.convasp.out',"w"),
            stderr = subprocess.PIPE
    ).communicate(stdin)
    if p[1] != "":
        print "\n" + p[1]
        print "See "+jobname+".convasp.out for more details"
        exit()
    

# === GULP ======================================================================

def prepareGulpInput(cifFile, gulpFile, jobName, verbose = False, auto_fractions = False): 
    """
    Example file:
    mgh2
    cell
    4.5168 4.5168 3.0205 90.0 90.0 90.0
    frac
    Mg 0.0 0.0 0.0
    H 0.306 0.306 0.0
    space
    136
    output xyz mgh2
    """
    from CifFile import CifFile
    
    if not os.path.exists(cifFile):
        raise IOError("CIF file '%s' was not found!" % (cifFile))
    
    cf = CifFile(cifFile)
    if verbose:
        print "------------------------------------------------------------------"
    if len(cf) != 1:
        raise StandardError("The cif file contains %i data blocks, while one was expected")
        # A cif file can contain several "datablocks" that each start
        # with "data_".
    
    if verbose:
        print "Reading data block '%s'..." % (cf.keys()[0])
    cb = cf[cf.keys()[0]]                               # open the first block
    AA = float(re.match('([0-9.e]*)',cb['_cell_length_a']).group(0))
    BB = float(re.match('([0-9.e]*)',cb['_cell_length_b']).group(0))
    CC = float(re.match('([0-9.e]*)',cb['_cell_length_c']).group(0))
    alpha = float(cb['_cell_angle_alpha'])
    beta = float(cb['_cell_angle_beta'])
    gamma = float(cb['_cell_angle_gamma'])    
    
    # Spacegroup number (1-230)
    # '_symmetry_Int_Tables_number' has been superseded by '_space_group_IT_number'
    
    if '_space_group_IT_number' in cb.keys():
        SG = int(cb['_space_group_IT_number'])
    elif '_symmetry_Int_Tables_number' in cb.keys():
        SG = int(cb['_symmetry_Int_Tables_number'])
    else:
        print "WARNING: No space group specified. Assuming P1."
        SG = 1
           

    # CCTBX:
    #unit_cell = uctbx.unit_cell([AA,BB,CC,alpha,beta,gamma])
    #space_group_info = sgtbx.space_group_info(symbol=cb['_symmetry_space_group_name_H-M'])
    #crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,space_group_info=space_group_info)    
    #print "CIF file read successfully:"
    #crystal_symmetry.show_summary()   

    if verbose:
        print "  Space group:",SG
        print "  a=%s, b=%s, c=%s, alpha=%s, beta=%s, gamma=%s" % (AA,BB,CC,alpha,beta,gamma)
    
    atomTypes = []
    atoms = ''
    fracOccFound = False
    firstAtom = True
    atoms = ""

    # The coordinates of the atom (_atom_site_fract_x/y/z) may have 
    # a last digit in parenthesis, like "-0.6636(7)". Therefore we
    # extract the part consisting of only digits and a decimal separator:
    coordsMatch = re.compile('[0-9.e-]*');

    for atom in cb.GetLoop('_atom_site_label'):
        atomKeys = dir(atom)
        if '_atom_site_type_symbol' in atomKeys:
            m = re.match('[a-z]*',atom._atom_site_type_symbol,re.I)
            atomType = m.group(0)
        else:
            m = re.match('[a-z]*',atom._atom_site_label,re.I)
            atomType = m.group(0)

        if '_atom_site_occupancy' in atomKeys:
            occ = float(coordsMatch.match(atom._atom_site_occupancy).group())
            if not occ == 1.0:
                if not fracOccFound: 
                    print " "
                print "  WARNING: Fractional occupancy (" + str(occ) +") " \
                    + "found for atom of type " + atomType + "."                
                fracOccFound = True
        else:
            occ = 1.0
        
        # Some crystal structures obtained by neutron diffraction use D for H:
        if atomType == 'D':
            atomType = 'H'
        
        if '_atom_site_symmetry_multiplicity' in atomKeys and '_atom_site_Wyckoff_symbol' in atomKeys:
            atomTypes.append(atomType+' at '+atom._atom_site_symmetry_multiplicity+atom._atom_site_Wyckoff_symbol)
        else:
            atomTypes.append(atomType)

        atomX = coordsMatch.match(atom._atom_site_fract_x).group()
        atomY = coordsMatch.match(atom._atom_site_fract_y).group()
        atomZ = coordsMatch.match(atom._atom_site_fract_z).group()
        
        atomPos = [atomX, atomY, atomZ]
        for i in range(3):
            pp = atomPos[i].split(".")
            if len(pp) is 2:
                decimals = pp[1]
                if len(decimals) > 3 and len(decimals) < 6 and decimals[0] == decimals[1] and decimals[0] == decimals[2] and decimals[-1] != "0":
                    if auto_fractions:
                        oldPos = atomPos[i]
                        atomPos[i] = "%.6f" % (float(eval('1.*'+float2fraction(atomPos[i]))))
                        print "  Notice: Converted %s into %s" %(oldPos,atomPos[i])
                    else:
                        print "\n"\
                            + "  ! Warning: The coordinate "+atomPos[i]+" looks similar to the fraction %s, but\n" % float2fraction(atomPos[i]) \
                            + "  !   has insufficient decimals to be recognized as so by GULP. If you want\n" \
                            + "  !   this coordinate to be recognized as a special high-symmetry position,\n" \
                            + "  !   you need to specify at least six digits. If you run cif2vasp with the \n" \
                            + "  !   -f switch, cif2vasp will try to add the necessary decimals automaticly."		
        
        atoms += "%s %s %s %s %f %f\n" % (atomType, atomPos[0], atomPos[1], atomPos[2], 0.0, occ)
        firstAtom = False

    if fracOccFound: 
        print " "
        print "ERROR: Fractional occupancies are not currently supported.\n"
        exit()
    
    if verbose:
        print "  Atom types: " + ', '.join(atomTypes)
    
    gulpFile = open(gulpFile,'w')       #Create and write the GULP  
    gulpFile.writelines([jobName+'\n',
        'cell\n',
        '%s %s %s %s %s %s\n' % (AA,BB,CC,alpha,beta,gamma),
        'frac\n',
        atoms,
        'space\n',
        str(SG)+'\n',
        'output xyz '+jobName+'\n'
    ])

# A probably not very robust function to convert a float like "0.333" to a fraction "1/3"
def float2fraction(f):
    f = float(f)
    num=1.                  # Start with 1 as numerator
    den=num/f               # Find denominator
    r=den%1
    if abs(r) > 1e-2:       # If denominator is decimal
        fac=1./r
        num *= fac          # Scale numerator..
        den *= fac          # .. and denominator
    
    return "%d/%d" % (round(num),round(den))


def runGulp(jobname, verbose = False):
    if verbose:
        print "Prepare primitive cell in xyz format using GULP... "
    p = subprocess.Popen("gulp", 
        stdin = file(jobname+'.gulp.in'), 
        stdout = file(jobname+'.gulp.out',"w"),
        stderr = subprocess.PIPE
    ).communicate() # communicate() returns a tuple (stdoutdata, stderrdata)
    if p[1] != "":
        print "\n" + p[1]
        print "See "+jobname+".gulp.out for more details"
        exit()

    if verbose:
        os.system("grep 'Crystal family' '%s'" % (jobname+'.gulp.out'))
        os.system("grep 'Space group' '%s'" % (jobname+'.gulp.out'))

if __name__ == "__main__":

    # Parse cmd line args
    parser = OptionParser( usage = "usage: %prog [options] filename.cif" )
    pdf_file = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
    parser.add_option('-v', '--verbose', action='store_true', dest = 'verbose', default = False, help = 'increase verbosity')
    parser.add_option('-f', action='store_true', dest = 'auto_fractions', default = False, help = 'try to identify fractional numbers (special positions) from coordinates with less than six decimals.')
    parser.add_option('-e', '--engine', dest = 'engine', default = 'ase', help = 'the backend to use for the conversion (gulp, ase or cctbx)')
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print "No filename given. Run %s -h for help" % (parser.get_prog_name())
        sys.exit(1)
    
    filename = args[0]
    if filename[-4:] == '.xyz':
        jobname = filename[0:-4]
        #cif2vaspUsingCCTBX(jobname)
        xyz2vaspUsingGULP(jobname, verbose = options.verbose, auto_fractions = options.auto_fractions)
    elif filename[-4:] == '.cif':
        jobname = filename[0:-4]
        if options.engine == 'cctbx':
            cif2vaspUsingCCTBX(jobname)
        elif options.engine == 'gulp':
            cif2vaspUsingGULP(jobname, verbose = options.verbose, auto_fractions = options.auto_fractions)
        elif options.engine == 'ase':
            cif2vaspUsingASE(jobname)
        else:
            print "Error: unknown engine specified."
    else:
        print "The input file must have the file-ending '.cif'"
        usage()
        sys.exit(2)
    
