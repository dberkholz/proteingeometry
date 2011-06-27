#!/usr/bin/env python

import sys
import optparse
from math import sqrt, cos

from Bio.PDB import *
import conformation_dependent_geometry.angles as angles

from shelxl_EandH_sigma import EandH_main_chain, EandH_side_chain

DFIX_format = 'DFIX_%s %.3f %.3f %s %s'
DANG_format = 'DANG_%s %.3f %.3f %s %s'
CHIV_format = 'CHIV_%s %.3f %s'

AminoAcids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'KCX',
              'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def calc_distance(atom1, atom2):
    """
    Calculates distance between atoms because this is not built into BIOPython

    scribed from http://www.scribd.com/doc/9816032/BioPython-for-Bioinfo
    """
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)

def angle_to_d13(angle, d12, d23):
    """ Given two bond lengths and the angle between them, what is the 
    distance between the end-points? """
    return sqrt(d12**2 + d23**2 - 2.0*d12*d23*cos(3.141592653/180*angle))

def sigma_of_d13(angle, sigma_angle, d12, d23):
    """ Assuming all uncertainty is in the value of the angle, what is the uncertainty of the 1-3 distance? """
    d13 = angle_to_d13(angle, d12, d23)
    try:
        sigma_d13 = (1.0/2.0)*(3.141592653/180.0)*sigma_angle/d13*sqrt(-(d12-d13-d23)*(d12+d13-d23)*(d12-d13+d23)*(d12+d13+d23))
    except ValueError:
        print 'We came up with a negative number in the square root in sigma_of_d13!'
        print angle, sigma_angle, d12, d23, d13
        sys.exit()
    return sigma_d13

def print_EandH_restraints(Residue, back_link, forward_link):
    """ Print the Engh & Huber restraints for Residue.  The target values 
    for the restraints that span parameterized as either PRO, GLY, or "general".
    The rest of the restraints come from the big tables in shelxl_EandH_sigma.
    The logical variables back_link and forward_link are true if the peptide
    chain continues in those directions. """
    residue_name=Residue['i'].get_id()[1]
    print 'REM Restraints for ', Residue['i'].get_resname()
    if back_link:
        if Residue['i'].get_resname() == 'PRO':
            print 'DFIX_{0} 1.341 0.016 C_- N'.format(residue_name)
            print 'DANG_{0} 2.463 0.059 C_- CA'.format(residue_name)
        elif Residue['i'].get_resname() == 'GLY':
            print 'DFIX_{0} 1.329 0.014 C_- N'.format(residue_name)
            print 'DANG_{0} 2.416 0.020 C_- CA'.format(residue_name)
        else:
            print 'DFIX_{0} 1.329 0.014 C_- N'.format(residue_name)
            print 'DANG_{0} 2.435 0.021 C_- CA'.format(residue_name)
    try:
        for restraint in EandH_main_chain[Residue['i'].get_resname()]:
            print restraint.format(residue_name)
        for restraint in EandH_side_chain[Residue['i'].get_resname()]:
            print restraint.format(residue_name)
    except:
        print 'Not implemented yet', Residue['i'].get_resname()
        exit
    if forward_link:
        # E&H 2001 is confusing but I think they say that these are the same for general, PRO, and GLY
        print 'DANG_{0} 2.425 0.026 CA N_+'.format(residue_name)
        print 'DANG_{0} 2.250 0.017 O N_+'.format(residue_name)
        print 'FLAT_{0} 0.5 O CA N_+ C CA_+'.format(residue_name)

def print_CDL_restraints(Residue, atoms, back_link, forward_link):
    """ Print the CDL restraints for Residue.  "atoms" contains the atoms for
    this and the adjoining residues.  The target values are parameterized by
    residue class and the phi/psi angles and extracted from the library using
    Donnie's conformation_dependent_geometry.angles.get_geometry.
    The logical variables back_link and forward_link are true if the peptide
    chain continues in those directions. """
    residue_name=Residue['i'].get_id()[1]
    phi = 180/3.141592653*calc_dihedral(atoms['C-'].get_vector(), atoms['N'].get_vector(), atoms['CA'].get_vector(), atoms['C'].get_vector())
    psi = 180/3.141592653*calc_dihedral(atoms['N'].get_vector(), atoms['CA'].get_vector(), atoms['C'].get_vector(), atoms['N+'].get_vector())
    print "REM  ", Residue['i'].get_resname(), Residue['i+1'].get_resname(), 'Phi = %.0f'%(phi), 'Psi = %.0f'%(psi)

    fields, geometry = angles.get_geometry(Residue['i'].get_resname(), Residue['i+1'].get_resname(), phi, psi)
    if back_link:
        print DFIX_format%(residue_name, geometry.L1Avg, geometry.L1Dev, 'C_-', 'N')
        print DANG_format%(residue_name, angle_to_d13(geometry.a1Avg,                 geometry.L1Avg, geometry.L2Avg), \
                                         sigma_of_d13(geometry.a1Avg, geometry.a1Dev, geometry.L1Avg, geometry.L2Avg), 'C_-', 'CA')

    print DFIX_format%(residue_name, geometry.L2Avg, geometry.L2Dev, 'N', 'CA')
    print DFIX_format%(residue_name, geometry.L4Avg, geometry.L4Dev, 'CA', 'C')
    print DFIX_format%(residue_name, geometry.L5Avg, geometry.L5Dev, 'C', 'O')

    print DANG_format%(residue_name, angle_to_d13(geometry.a3Avg,                 geometry.L2Avg, geometry.L4Avg), \
                                     sigma_of_d13(geometry.a3Avg, geometry.a3Dev, geometry.L2Avg, geometry.L4Avg), 'N', 'C')
    print DANG_format%(residue_name, angle_to_d13(geometry.a5Avg,                 geometry.L4Avg, geometry.L5Avg), \
                                     sigma_of_d13(geometry.a5Avg, geometry.a5Dev, geometry.L4Avg, geometry.L5Avg), 'CA', 'O')
    print 'CHIV_{0} C'.format(residue_name)
    if Residue['i'].get_resname() != 'GLY': 
        print DFIX_format%(residue_name, geometry.L3Avg, geometry.L3Dev, 'CA', 'CB')
        print DANG_format%(residue_name, angle_to_d13(geometry.a2Avg,                 geometry.L2Avg, geometry.L3Avg), \
                                         sigma_of_d13(geometry.a2Avg, geometry.a2Dev, geometry.L2Avg, geometry.L3Avg), 'N', 'CB')
        print DANG_format%(residue_name, angle_to_d13(geometry.a4Avg,                 geometry.L4Avg, geometry.L3Avg), \
                                         sigma_of_d13(geometry.a4Avg, geometry.a4Dev, geometry.L4Avg, geometry.L3Avg), 'C', 'CB')
        print CHIV_format%(residue_name, geometry.L2Avg*geometry.L3Avg*geometry.L4Avg* \
                    sqrt(1.0 - cos(3.141592653/180*geometry.a2Avg)**2 - \
                               cos(3.141592653/180*geometry.a3Avg)**2 - \
                               cos(3.141592653/180*geometry.a4Avg)**2 + 
                           2.0*cos(3.141592653/180*geometry.a2Avg)* \
                               cos(3.141592653/180*geometry.a3Avg)* \
                               cos(3.141592653/180*geometry.a4Avg)), 'CA')

    if Residue['i'] not in ('GLY', 'ALA'):
        for restraint in EandH_side_chain[Residue['i'].get_resname()]:
            print restraint.format(residue_name)

    if forward_link:
        # Length of C-N_+ is a problem because it doesn't come with this residue in the CDL.
        # We'll use E&H 1991 for this length.  Proline is a special case.
        if Residue['i+1'].get_resname() == 'PRO':
            CNlen = 1.341
        else: 
            CNlen = 1.329
        print DANG_format%(residue_name, angle_to_d13(geometry.a6Avg,                 geometry.L4Avg, CNlen), \
                                         sigma_of_d13(geometry.a6Avg, geometry.a6Dev, geometry.L4Avg, CNlen), 'CA', 'N_+')
        print DANG_format%(residue_name, angle_to_d13(geometry.a7Avg,                 geometry.L5Avg, CNlen), \
                                         sigma_of_d13(geometry.a7Avg, geometry.a7Dev, geometry.L5Avg, CNlen), 'O', 'N_+')
        print 'FLAT_{0} 0.5 O CA N_+ C CA_+'.format(residue_name)

def main():
    """This is the main coordinating routine for this package.
    It reads the PDB file, loops over the residues of each chain 
    and writes out the CDL restraints in Shelxl format."""

    args_parser = optparse.OptionParser(version='%prog ' + '0.1')
    args_parser.set_usage('usage: cdl-shelxl.py <PDB filename>')
    optlist, args = args_parser.parse_args()

    if len(args) != 1:
        args_parser.error('Not the right number of arguments')
    else:
        PDBfilename = args[0]

    parser=PDBParser()
    s=parser.get_structure('PDB', PDBfilename)

    chain=s[0].get_list()
    if len(chain) != 1:
        print 'You do realize that shelxl cannot handle models with more than one chain!'
        sys.exit()
    print 'REM Working with chain ', s[0].id, '/', chain[0].id

    Residue = {}
    #Residue['i+1'] = None
    State = 'Empty'
    #print chain[0].get_iterator().next()
    chain_iterator = chain[0].get_iterator()

    while State != 'Finished':
       if State == 'Empty':
          Residue['i-1'] = None
          Residue['i']   = None
          try:
             Residue['i+1']=chain_iterator.next()
          except StopIteration:
             State = 'Finished'
             break
          #print Residue['i+1'].get_resname()
          if Residue['i+1'].get_resname() in AminoAcids:         
             State = 'Filling'
          #else:
             #print 'State: ', State, '\n    Skipping residue', Residue['i+1']

       elif State == 'Filling':
          #Residue['i-1']=None
          Residue['i'] = Residue['i+1']
          try:
             Residue['i+1']=chain_iterator.next()
          except StopIteration:
             print '    exception!'
             Residue['i+1'] = None
          #print 'State: ', State, '\n    Working on residue', Residue['i'].get_id(), ' as E&H(3).'
          print_EandH_restraints(Residue, False, Residue['i+1']!=None)
          if Residue['i+1'] == None or Residue['i+1'].get_resname() not in AminoAcids:
             State = 'Empty'
             Residue['i+1'] = None
          else: State = 'Running'

       elif State == 'Running':
          Residue['i-1']=Residue['i']
          Residue['i']=Residue['i+1']
          try:
             Residue['i+1']=chain_iterator.next()
             if Residue['i+1'].get_resname() not in AminoAcids:
                 Residue['i+1']=None
                 State = 'Empty'
          except StopIteration:
             Residue['i+1']=None
             State = 'Empty'
          atoms = {}
          for atom in Residue['i-1'].get_unpacked_list():
             if atom.name == 'C':
                atoms['C-'] = atom
          for atom in Residue['i'].get_unpacked_list():
             atoms[atom.name] = atom
          if Residue['i+1'] != None:
             for atom in Residue['i+1'].get_unpacked_list():
                if atom.name == 'N':
                   atoms['N+'] = atom
          back_link    = (Residue['i-1'] != None and 
                          Residue['i-1'].get_resname() in AminoAcids and
                          calc_distance(atoms['C-'].get_vector(), atoms['N' ].get_vector()) < 2.5)
          forward_link = (Residue['i+1'] != None and 
                          Residue['i+1'].get_resname() in AminoAcids and
                          calc_distance(atoms['C' ].get_vector(), atoms['N+'].get_vector()) < 2.5)
          # Can we be in state=running if back_link is false?  Yes, because we could have a break or a non-AA residue.
          if back_link and forward_link:
             #print 'State: ', State, '\n    Working on CDL for residue', Residue['i'].get_id()
             print_CDL_restraints(Residue, atoms, back_link, forward_link)
          else:
             #print 'State: ', State, '\n    Working on residue', Residue['i'].get_id(), ' as E&H(2).'
             print_EandH_restraints(Residue, back_link, forward_link)

# I expect this code will always be run as a program, but just in case.
if __name__ == "__main__":
   sys.exit(main())

