#!/usr/bin/env python

import sys
import optparse
from math import sqrt, cos

from Bio.PDB import *
import conformation_dependent_geometry.angles as angles

from conformation_dependent_geometry.shelxl_EandH_sigma \
     import EandH_main_chain, EandH_side_chain

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

def EandH_restraints(Residue, back_link, forward_link):
    """ Create a list of strings with the Engh & Huber restraints for Residue.
    The target values for the restraints that span parameterized as either
    PRO, GLY, or "general".
    The rest of the restraints come from the big tables in shelxl_EandH_sigma.
    The logical variables back_link and forward_link are true if the peptide
    chain continues in those directions. """
    restraint_list=[]
    residue_name=Residue['i'].get_id()[1]
    restraint_list.append('REM Restraints for '+Residue['i'].get_resname())
    if back_link:
        if Residue['i'].get_resname() == 'PRO':
            restraint_list.append('DFIX_{0} 1.341 0.016 C_- N'.format(residue_name))
            restraint_list.append('DANG_{0} 2.463 0.059 C_- CA'.format(residue_name))
        elif Residue['i'].get_resname() == 'GLY':
            restraint_list.append('DFIX_{0} 1.329 0.014 C_- N'.format(residue_name))
            restraint_list.append('DANG_{0} 2.416 0.020 C_- CA'.format(residue_name))
        else:
            restraint_list.append('DFIX_{0} 1.329 0.014 C_- N'.format(residue_name))
            restraint_list.append('DANG_{0} 2.435 0.021 C_- CA'.format(residue_name))
    try:
        for restraint in EandH_main_chain[Residue['i'].get_resname()]:
            restraint_list.append(restraint.format(residue_name))
        for restraint in EandH_side_chain[Residue['i'].get_resname()]:
            restraint_list.append(restraint.format(residue_name))
    except:
        restraint_list.append('Not implemented yet'+Residue['i'].get_resname())
        exit
    if forward_link:
        # E&H 2001 is confusing but I think they say that these are the same for general, PRO, and GLY
        restraint_list.append('DANG_{0} 2.425 0.026 CA N_+'.format(residue_name))
        restraint_list.append('DANG_{0} 2.250 0.017 O N_+'.format(residue_name))
        restraint_list.append('FLAT_{0} 0.5 O CA N_+ C CA_+'.format(residue_name))
    return restraint_list

def CDL_restraints(Residue, atoms, back_link, forward_link):
    """ Print the CDL restraints for Residue.  "atoms" contains the atoms for
    this and the adjoining residues.  The target values are parameterized by
    residue class and the phi/psi angles and extracted from the library using
    Donnie's conformation_dependent_geometry.angles.get_geometry.
    The logical variables back_link and forward_link are true if the peptide
    chain continues in those directions. """
    restraint_list=[]
    residue_name=Residue['i'].get_id()[1]
    phi = 180/3.141592653*calc_dihedral(atoms['C-'].get_vector(), atoms['N'].get_vector(), atoms['CA'].get_vector(), atoms['C'].get_vector()) + 180
    psi = 180/3.141592653*calc_dihedral(atoms['N'].get_vector(), atoms['CA'].get_vector(), atoms['C'].get_vector(), atoms['N+'].get_vector()) + 180
    restraint_list.append("REM  "+Residue['i'].get_resname()+" "+Residue['i+1'].get_resname()+' Phi = %.0f'%(phi)+' Psi = %.0f'%(psi))
    fields, geometry = geom[Residue['i'].get_resname(), Residue['i+1'].get_resname(), phi, psi]
    if back_link:
        restraint_list.append(DFIX_format%(residue_name, getattr(geometry, 'C(-1)-NAvg(i)'), getattr(geometry, 'C(-1)-NDev(i)'), 'C_-', 'N'))
        restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'C(-1)-N-CAAvg(i)'),                 getattr(geometry, 'C(-1)-NAvg(i)'), getattr(geometry, 'N-CAAvg(i)')), \
                                         sigma_of_d13(getattr(geometry, 'C(-1)-N-CAAvg(i)'), getattr(geometry, 'C(-1)-N-CADev(i)'), getattr(geometry, 'C(-1)-NAvg(i)'), getattr(geometry, 'N-CAAvg(i)')), 'C_-', 'CA'))

    restraint_list.append(DFIX_format%(residue_name, getattr(geometry, 'N-CAAvg(i)'), getattr(geometry, 'N-CADev(i)'), 'N', 'CA'))
    restraint_list.append(DFIX_format%(residue_name, getattr(geometry, 'CA-CAvg(i)'), getattr(geometry, 'CA-CDev(i)'), 'CA', 'C'))
    restraint_list.append(DFIX_format%(residue_name, getattr(geometry, 'C-OAvg(i)'), getattr(geometry, 'C-ODev(i)'), 'C', 'O'))

    restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'N-CA-CAvg(i)'),                 getattr(geometry, 'N-CAAvg(i)'), getattr(geometry, 'CA-CAvg(i)')), \
                                     sigma_of_d13(getattr(geometry, 'N-CA-CAvg(i)'), getattr(geometry, 'N-CA-CDev(i)'), getattr(geometry, 'N-CAAvg(i)'), getattr(geometry, 'CA-CAvg(i)')), 'N', 'C'))
    restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'CA-C-OAvg(i)'),                 getattr(geometry, 'CA-CAvg(i)'), getattr(geometry, 'C-OAvg(i)')), \
                                     sigma_of_d13(getattr(geometry, 'CA-C-OAvg(i)'), getattr(geometry, 'CA-C-ODev(i)'), getattr(geometry, 'CA-CAvg(i)'), getattr(geometry, 'C-OAvg(i)')), 'CA', 'O'))
    restraint_list.append('CHIV_{0} C'.format(residue_name))
    if Residue['i'].get_resname() != 'GLY': 
        restraint_list.append(DFIX_format%(residue_name, getattr(geometry, 'CA-CBAvg(i)'), getattr(geometry, 'CA-CBDev(i)'), 'CA', 'CB'))
        restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'N-CA-CBAvg(i)'),                 getattr(geometry, 'N-CAAvg(i)'), getattr(geometry, 'CA-CBAvg(i)')), \
                                         sigma_of_d13(getattr(geometry, 'N-CA-CBAvg(i)'), getattr(geometry, 'N-CA-CBDev(i)'), getattr(geometry, 'N-CAAvg(i)'), getattr(geometry, 'CA-CBAvg(i)')), 'N', 'CB'))
        restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'CB-CA-CAvg(i)'),                 getattr(geometry, 'CA-CAvg(i)'), getattr(geometry, 'CA-CBAvg(i)')), \
                                         sigma_of_d13(getattr(geometry, 'CB-CA-CAvg(i)'), getattr(geometry, 'CB-CA-CDev(i)'), getattr(geometry, 'CA-CAvg(i)'), getattr(geometry, 'CA-CBAvg(i)')), 'C', 'CB'))
        restraint_list.append(CHIV_format%(residue_name, getattr(geometry, 'N-CAAvg(i)')*getattr(geometry, 'CA-CBAvg(i)')*getattr(geometry, 'CA-CAvg(i)')* \
                    sqrt(1.0 - cos(3.141592653/180*getattr(geometry, 'N-CA-CBAvg(i)'))**2 - \
                               cos(3.141592653/180*getattr(geometry, 'N-CA-CAvg(i)'))**2 - \
                               cos(3.141592653/180*getattr(geometry, 'CB-CA-CAvg(i)'))**2 + 
                           2.0*cos(3.141592653/180*getattr(geometry, 'N-CA-CBAvg(i)'))* \
                               cos(3.141592653/180*getattr(geometry, 'N-CA-CAvg(i)'))* \
                               cos(3.141592653/180*getattr(geometry, 'CB-CA-CAvg(i)'))), 'CA'))

    if Residue['i'] not in ('GLY', 'ALA'):
        for restraint in EandH_side_chain[Residue['i'].get_resname()]:
            restraint_list.append(restraint.format(residue_name))

    if forward_link:
        # Length of C-N_+ is a problem because it doesn't come with this residue in the CDL.
        # We'll use E&H 1991 for this length.  Proline is a special case.
        if Residue['i+1'].get_resname() == 'PRO':
            CNlen = 1.341
        else: 
            CNlen = 1.329
        restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'CA-C-N(+1)Avg(i)'),                 getattr(geometry, 'CA-CAvg(i)'), CNlen), \
                                         sigma_of_d13(getattr(geometry, 'CA-C-N(+1)Avg(i)'), getattr(geometry, 'CA-C-N(+1)Dev(i)'), getattr(geometry, 'CA-CAvg(i)'), CNlen), 'CA', 'N_+'))
        restraint_list.append(DANG_format%(residue_name, angle_to_d13(getattr(geometry, 'O-C-N(+1)Avg(i)'),                 getattr(geometry, 'C-OAvg(i)'), CNlen), \
                                         sigma_of_d13(getattr(geometry, 'O-C-N(+1)Avg(i)'), getattr(geometry, 'O-C-N(+1)Dev(i)'), getattr(geometry, 'C-OAvg(i)'), CNlen), 'O', 'N_+'))
        restraint_list.append('FLAT_{0} 0.5 O CA N_+ C CA_+'.format(residue_name))
    return restraint_list

def process(PDBfilename, format):
    """This is the main coordinating routine for this package.
    It reads the PDB file, loops over the residues of each chain 
    and writes out the CDL restraints in Shelxl format."""

    if format == 'text':
        terminus = ''
    elif format == 'html':
        terminus = '<br>'        
    else:
        print 'Second parameter must be either \'test\' or \'html\'.'
        return 1

    parser=PDBParser()
    s=parser.get_structure('PDB', PDBfilename)

    chain=s[0].get_list()
    if len(chain) != 1:
        print 'You do realize that <i>shelxl</i> cannot handle models with more than one chain!', terminus
        return 1
    if chain[0].id != ' ':
        print 'Your residues are labeled as chain ', chain[0].id, '.  shelxl requires the chain id to be blank.', terminus
        return 1

    if format == 'html':
        print "Use copy-and-paste to move these lines into your .ins file, replacing the Engh & Huber definitions.", terminus, terminus

    global geom
    geom = angles.setup()

    Residue = {}
    #Residue['i+1'] = None
    State = 'Empty'
    #print chain[0].get_iterator().next(), terminus
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
          #print Residue['i+1'].get_resname(), terminus
          if Residue['i+1'].get_resname() in AminoAcids:         
             State = 'Filling'
          #else:
             #print 'State: ', State, '\n    Skipping residue', Residue['i+1'], terminus

       elif State == 'Filling':
          #Residue['i-1']=None
          Residue['i'] = Residue['i+1']
          try:
             Residue['i+1']=chain_iterator.next()
          except StopIteration:
             print '    exception!', terminus
             Residue['i+1'] = None
          #print 'State: ', State, '\n    Working on residue', Residue['i'].get_id(), ' as E&H(3).', terminus
          for line in EandH_restraints(Residue, False, Residue['i+1']!=None):
              print line, terminus
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
             #print 'State: ', State, '\n    Working on CDL for residue', Residue['i'].get_id(), terminus
             for line in CDL_restraints(Residue, atoms, back_link, forward_link):
                 print line, terminus
          else:
             #print 'State: ', State, '\n    Working on residue', Residue['i'].get_id(), ' as E&H(2).', terminus
             for line in EandH_restraints(Residue, back_link, forward_link):
                 print line, terminus
    return 0

# Currently this code is run directly when producing a text file but "process" is called directly to produce html.
if __name__ == "__main__":

    args_parser = optparse.OptionParser(version='%prog ' + '0.1')
    args_parser.set_usage('usage: cdl-shelxl.py <PDB filename>')
    optlist, args = args_parser.parse_args()

    if len(args) != 1:
        args_parser.error('Not the right number of arguments')
    else:
        PDBfilename = args[0]

    sys.exit(process(PDBfilename, 'text'))

