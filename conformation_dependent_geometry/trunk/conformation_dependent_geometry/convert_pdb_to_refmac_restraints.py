#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright © 2007-2008 Oregon State University
# Copyright © 2010 Mayo Clinic College of Medicine
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
# Authors:
#     Donnie Berkholz <donnie.berkholz@gmail.com>

import math, sys
# Do it this way so __init__.py runs and initializes angles.optlist. Otherwise
# we get tracebacks.
import conformation_dependent_geometry.angles as angles

try:
    import mmLib.Structure, mmLib.FileIO, mmLib.AtomMath, mmLib.ConsoleOutput
except:
    print "Failed to import mmLib"
    print "Install from http://pymmlib.sourceforge.net/"
    sys.exit(1)

# Program version                                                              
version = '0.1'

def print_restraints(dbdict, struct):
    # Basic design:
    # - Read PDB
    # - Iterate over all residues
    #   - For each residue:
    #     - Find phi/psi
    #     - Check for atoms in our angles and look up ideal values

    # First, build the properties. Since we need lookahead to find
    # conformation-dependent restraints, we'll iterate a second time to
    # simplify the code
    for r in struct.iter_amino_acids():
        res_seq, icode = mmLib.Structure.fragment_id_split(r.fragment_id)

        r.phi = r.calc_torsion_phi()
        r.psi = r.calc_torsion_psi()
        if not r.phi:
            r.phi = math.radians(999)
        if not r.psi:
            r.psi = math.radians(999)

        r.props = {
            'name': r.res_name,
            'res_num': res_seq,
            'icode': icode,
            'chain': r.chain_id,
            'phi': math.degrees(r.phi),
            'psi': math.degrees(r.psi),
        }

    # Get the conformation-dependent values
    for r in struct.iter_amino_acids():
        if r.props['phi'] == 999 or r.props['psi'] == 999:
            continue

        next_res = r.get_offset_residue(1)
        prev_res = r.get_offset_residue(-1)
        prev_res_seq, prev_icode = mmLib.Structure.fragment_id_split(
            prev_res.fragment_id)
        res_seq, icode = mmLib.Structure.fragment_id_split(
            r.fragment_id)
        next_res_seq, next_icode = mmLib.Structure.fragment_id_split(
            next_res.fragment_id)

        # Make sure the adjacent residues really are sequential, otherwise we
        # can get weird angles.  Also account for insertion codes, where
        # res_seq is equal but icode changes
        if next_res:
            if (int(next_res_seq) - int(res_seq) != 1):
                if (int(next_res_seq) - int(res_seq) != 0) \
                       and next_icode and icode \
                       and next_icode != icode:
                    pass
                else:
                    continue
        if prev_res:
            if (int(prev_res_seq) - int(res_seq) != -1):
                if (int(prev_res_seq) - int(res_seq) != 0) \
                       and prev_icode and icode \
                       and prev_icode != icode:
                    pass
                else:
                    continue

        fields, geometry = angles.get_geometry(
                             dbdict, 
                             r.props['name'], 
                             next_res.props['name'], 
                             r.props['phi'], 
                             r.props['psi'])

        for angle_definition, angle_name in zip(
            angles.angle_atoms,
            angles.angle_names):

            atoms = []
            for residue_index, atom_name in angle_definition:
                # check for residue existence
                # check for atom existence
                if residue_index == 0:
                    atom = prev_res.get_atom(atom_name)
                elif residue_index == 1:
                    atom = r.get_atom(atom_name)
                elif residue_index == 2:
                    atom = next_res.get_atom(atom_name)
                atoms.append(atom)
            if len(atoms) != 3 or None in atoms:
                # the angle is missing atoms, so don't restrain it
                continue

            param_name = angle_name

            pgd_param_average_name = \
                angles.get_database_attribute_average_name(param_name)
            pgd_param_deviation_name = \
                angles.get_database_attribute_deviation_name(param_name)
            # provide a geometry object
            new_angle_ideal = getattr(geometry, pgd_param_average_name)
            new_deviation_ideal = getattr(geometry, pgd_param_deviation_name)

            print 'external angle',
            for i in range(0, 3):
                if i == 0:
                    atom_position = 'first'
                elif i == 1:
                    atom_position = 'second'
                elif i == 2:
                    atom_position = 'third'
                print atom_position, 'atom', atoms[i].name,
                if atoms[i].alt_loc:
                    print 'alt', atoms[i].alt_loc,
                if atoms[i].chain_id:
                    print 'chain', atoms[i].chain_id,
                res_seq, icode = mmLib.Structure.fragment_id_split(
                    atoms[i].fragment_id)
                print 'residue', res_seq,
                if icode:
                    print 'insertion', icode,

            print 'value', new_angle_ideal,
            print 'sigma', new_deviation_ideal,
            print 'type 0'


def optparse_setup():
    import optparse

    usage = """usage: %prog [options] [<PDB>]"""

    parser = optparse.OptionParser(version='%prog ' + version)
    parser.set_usage(usage)
    parser.set_defaults(
        verbose=False,
        )
    parser.add_option( \
        '-v', '--verbose', \
            action='store_true', \
            dest='verbose', \
            help='Use verbose output; defaults to %default')
    return parser


def main(argv):
    global optlist
    global args
    parser = optparse_setup()
    optlist, args = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    dbdict = angles.create_all_databases(angles.databases)


    # Shut up mmLib before we load any structures
    mmLib.ConsoleOutput.disable()

    pdb = args[0]
    struct = mmLib.FileIO.LoadStructure(file=pdb)
    print_restraints(dbdict, struct)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
