#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright © 2007-2008 Oregon State University
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
#     Donnie Berkholz <berkhold@science.oregonstate.edu>

import math, sys
# Do it this way so __init__.py runs and initializes angles.optlist. Otherwise
# we get tracebacks.
import conformation_dependent_geometry.angles as angles

try:
    import mmLib.Structure, mmLib.FileIO, mmLib.AtomMath
except:
    print "Failed to import mmLib"
    print "Install from http://pymmlib.sourceforge.net/"
    sys.exit(1)

# Structure.AminoAcidResidue:
#
# calc_mainchain_bond_angle(self)
# Calculates main chain bond angles (N-CA-C, N-CA-CB, CB-CA-C, CA-C-O,
# CA-C-(next)N, C-(next residue)N-(next residue)CA) and returns the result as
# a 6-tuple in that order.
#
# calc_mainchain_bond_length(self)
# Calculates the main chain bond lengths: (N-CA, CA-C, C-O, CA-CB, CA-(next)N).
#
# calc_torsion_phi(self)
# Calculates the Phi torsion angle of the amino acid.
#
# calc_torsion_psi(self)
# Calculates the Psi torsion angle of the amino acid.
#
# calc_torsion_omega(self)
# Calculates the Omega torsion angle of the amino acid.

# Program version                                                              
version = '0.1'

class geom(object):
    def __init__(self, attr):
        setattr(self, attr, 0)
    def set(self, attr, value):
        setattr(self, attr, value)
    def get(self, attr):
        return getattr(self, attr)

class protein_geometry_database(geom):
    def __init__(self, attr):
        geom.__init__(self, attr)
        # Read in database files etc here
        self.dbdict = angles.create_all_databases(
                        angles.databases)

    def get(self, phi, psi, attr):
        # Get based on residue type, phi and psi
        #print phi, psi
        #print type(attr)

        pgdattr = angles.get_database_attribute_average_name(attr)

        if phi > 180 or psi > 180:
            residue = 'default'
            phi = 0
            psi = 0
        # print 'residue = ', residue
        fields, geometry = angles.get_geometry(
                             self.dbdict, 
                             self.res_name, 
                             self.next_res_name, 
                             phi, 
                             psi)
        #print geometry.__dict__
        value = getattr(geometry, pgdattr)
        obs = getattr(geometry, 'Observations')
        # print 'res, phi, psi, value, observations = %s, %.2f, %.2f, %.2f, %d' % (residue, phi, psi, value, obs)
        # Our defaults don't include omega for some reason,
        # instead they return -1
        if attr == 'ome' and value == -1:
            return 180
        else:
            return value

class CisPeptide(Exception):
    def __str__(self):
        return 'caught %s' % repr(self)

class BadDihedral(Exception):
    def __str__(self):
        return 'caught %s' % repr(self)

class NoEquivalentFragment(Exception):
    def __str__(self):
        return 'caught %s' % repr(self)

class NoMeasurement(Exception):
    def __str__(self):
        return 'caught %s' % repr(self)

class InvalidDeviation(Exception):
    def __str__(self):
        return 'caught %s' % repr(self)

class pdb_container(object):
    def __init__(self, name):
        self.name = name
#    def __repr__(self):
#        return self.name
    def __str__(self):
        return self.name

def valid_deviation(geomclass, deviation):
    """If a deviation is too large to be reasonable, ignore it"""
    max_angle_deviation = 40
    max_length_deviation = 0.5
    max_torsion_deviation = 60
    deviation = abs(deviation)

    if geomclass == 'angle':
        if deviation > max_angle_deviation:
            return False
    if geomclass == 'length':
        if deviation > max_length_deviation:
            return False
    if geomclass == 'torsion':
        if deviation > max_torsion_deviation:
            return False
    return True

def main(argv):
    global optlist
    global args
    parser = optparse_setup()
    optlist, args = parser.parse_args()

    if len(args) != 2:
        parser.error('incorrect number of arguments')

    meas = args[0]

    pdb = args[1]
    struct = mmLib.FileIO.LoadStructure(file=pdb)
    get_geometry(struct, meas)

    # To make validation easier
    if meas.startswith('a'):
        geomclass = 'angle'
    elif meas.startswith('l'):
        geomclass = 'length'
    elif meas == 'ome' \
            or meas == 'zeta':
        geomclass = 'torsion'

    if optlist.compare_pdbs:
        cpdbs = []
        cstructs = {}
        msd = {}
        for cpdb in optlist.compare_pdbs:
            cpdbs.append(pdb_container(name=cpdb))
        for cpdb in cpdbs:
            cstructs[cpdb] = mmLib.FileIO.LoadStructure(file=cpdb.name)
            get_geometry(cstructs[cpdb], meas)
            msd[cpdb] = 0
    
    if optlist.compare_eh:
        eh = geom(meas)
        # print eh.__dict__
        if meas == 'ome':
            # print "Measurement is ome"
            eh.set(meas, 180)
        eh.msd = 0

    if optlist.compare_cdecg:
        cdecg = protein_geometry_database(meas)
        cdecg.msd = 0

    N = 0
    for r in struct.iter_amino_acids():
        r_atom = r.get_atom('CA')

        if optlist.compare_pdbs:
            c_atoms = {}
            c_residues = {}
            try:
                for cpdb in cstructs.keys():
                    c_atoms[cpdb] = cstructs[cpdb].get_equivalent_atom(r_atom)
                    # c_residues are the compared residues
                    try:
                        c_residues[cpdb] = c_atoms[cpdb].get_fragment()
                    except:
                        raise NoEquivalentFragment
                    try:
                        c_residues[cpdb].props[meas]
                    except:
                        raise NoMeasurement
            except (NoEquivalentFragment, NoMeasurement), e:
                print e
                continue
        try:
            r.props[meas]
        except:
            continue

        if optlist.compare_cdecg:
            # We use this to look up conformation-dependent geometry
            cdecg.res_name = r.res_name
            # Look for residue i+1 for a proline, so we know if we're XPro
            res_seq, icode = mmLib.Structure.fragment_id_split(r.fragment_id)
            cdecg.next_res_numstr = str(int(res_seq) + 1)
            r_dict = r.get_chain().fragment_dict
            try:
                cdecg.next_res_name = r_dict[cdecg.next_res_numstr]
            except:
                cdecg.next_res_name = 'END'

        if optlist.compare_eh:
            if meas == 'a3':
                # print "Measurement is a3"
                if r.res_name == 'GLY':
                    # print r2
                    eh.set(meas, 112.5)
                elif r.res_name == 'PRO':
                    # print r2
                    eh.set(meas, 111.8)
                else:
                    eh.set(meas, 111.2)
        if meas == 'phi' or meas == 'psi' or meas == 'ome':
            #print "phi/psi/ome"
            if r.props[meas] > 180:
                continue
            if optlist.compare_pdbs:
                # if _any_ residues in set are bad, skip comparing this
                # residue altogether (XXX FIXME: is this a good idea?)
                try:
                    for cpdb in c_residues.keys():
                        if c_residues[cpdb].props[meas] > 180:
                            raise BadDihedral
                except BadDihedral, e:
                    print e
                    continue
        # Force omega into the -90 to +270 range so it's an easy comparison
        # The other option is doing math in radians with a pi modulus
        if meas == 'ome':
            if r.props[meas] < -90:
                r.props[meas] += 360
            if optlist.compare_pdbs:
                for cpdb in c_residues.keys():
                    if c_residues[cpdb].props[meas] < -90:
                        c_residues[cpdb].props[meas] += 360
            # ignore cis peptides
            try:
                if r.props[meas] < 90 and r.props[meas] > -90:
                    raise CisPeptide
                if optlist.compare_pdbs:
                    for cpdb in c_residues.keys():
                        if c_residues[cpdb].props[meas] < 90 and c_residues[cpdb].props[meas] > -90:
                            raise CisPeptide
            except CisPeptide, e:
                print e
                continue

        if optlist.compare_pdbs:
            dev = {}
            try:
                for cpdb in cstructs.keys():
                    dev[cpdb] = r.props[meas] - c_residues[cpdb].props[meas]
                    if not valid_deviation(geomclass, dev[cpdb]):
                        raise InvalidDeviation
                    msd[cpdb] += dev[cpdb]**2
            except InvalidDeviation, e:
                print e
                continue
        if optlist.compare_eh:
            eh.dev = r.props[meas] - eh.get(meas)
            if not valid_deviation(geomclass, eh.dev):
                continue
            eh.msd += eh.dev**2
        if optlist.compare_cdecg:
            cdecg_meas = cdecg.get(r.props['phi'], r.props['psi'], meas)
            cdecg.dev = r.props[meas] - cdecg_meas
            if not valid_deviation(geomclass, cdecg.dev):
                continue
            cdecg.msd += cdecg.dev**2
        N += 1
        if optlist.verbose:
            print '%s %s %s' % (r.props['name'], r.props['chain'], r.props['id']),
            print '%+.1f %+.1f' % (r.props['phi'], r.props['psi']),
            print '%.2f' % r.props[meas],
            if optlist.compare_pdbs:
                for cpdb in c_residues.keys():
                    print '%.2f' % c_residues[cpdb].props[meas],
            if optlist.compare_eh:
                print '%.2f' % eh.get(meas),
            if optlist.compare_cdecg:
                print '%.2f' % cdecg_meas,
            if optlist.compare_pdbs:
                for cpdb in cstructs.keys():
                    print '%+.2f %.2f' % (dev[cpdb], msd[cpdb]),
            if optlist.compare_eh:
                print '%+.2f %.2f' % (eh.dev, eh.msd),
            if optlist.compare_cdecg:
                print '%+.2f %.2f' % (cdecg.dev, cdecg.msd),
            print '%d' % N
    print
    print 'Using %d residues' % N
    if optlist.compare_pdbs:
        rmsd = {}
        for cpdb in cstructs.keys():
            rmsd[cpdb] = math.sqrt ( msd[cpdb] / N )
            print '%s for %s vs %s = %.2f' % (meas, pdb, cpdb, rmsd[cpdb])
    if optlist.compare_eh:
        eh.rmsd = math.sqrt ( eh.msd / N )
        print '%s for %s vs E&H = %.2f' % (meas, pdb, eh.rmsd)
    if optlist.compare_cdecg:
        cdecg.rmsd = math.sqrt ( cdecg.msd / N )
        print '%s for %s vs CD-ECG = %.2f' % (meas, pdb, cdecg.rmsd)

def get_geometry(struct, *geomtypes):
    # calculate bond angles
    angle=False
    # calculate bond lengths
    length=False
    # calculate torsions besides phi and psi
    torsion=False
    for geomtype in geomtypes:
        if geomtype.startswith('a'):
            angle=True
        elif geomtype.startswith('l'):
            length=True
        elif geomtype == 'ome' \
                or geomtype == 'zeta':
            torsion=True

    for r in struct.iter_amino_acids():
        next_res = r.get_offset_residue(1)
        prev_res = r.get_offset_residue(-1)

        # Make sure the adjacent residues really are sequential
        # Otherwise we can get weird angles
        if next_res:
            if (int(next_res.fragment_id) - int(r.fragment_id) != 1):
                continue
        if prev_res:
            if (int(prev_res.fragment_id) - int(r.fragment_id) != -1):
                continue

        if angle:
            r.a3, r.a2, r.a4, r.a5, a6, a1 = \
                r.calc_mainchain_bond_angle()
            if not prev_res:
                # n-terminus
                r.a1 = 0
                r.a6 = 0
                r.a7 = 0
        if length:
            r.l2, r.l4, r.l5, r.l3, r.l6 = \
                r.calc_mainchain_bond_length()
            r.l7 = 0
            if prev_res:
                r.l1 = prev_res.l6
            else:
                # n-terminus
                r.l1 = 0
        if next_res:
            if angle:
                next_res.a1 = a1
                next_res.a6 = a6
                aO = r.get_atom('O')
                aC = r.get_atom('C')
                naN = next_res.get_atom('N')
                next_res.a7 = mmLib.AtomMath.calc_angle(aO, aC, naN)
            if length:
                aN = r.get_atom('N')
                aCA = r.get_atom('CA')
                next_res.l2 = mmLib.AtomMath.calc_distance(aN, aCA)
                r.l7 = next_res.l2

        else:
            if length:
                # c-terminus
                r.l7 = 0
        r.phi = r.calc_torsion_phi()
        r.psi = r.calc_torsion_psi()
        if torsion:
            r.ome = r.calc_torsion_omega()

        # glycine
        if length:
            if not r.l3:
                r.l3 = 0
        if angle:
            if not r.a2:
                r.a2 = 0
                r.a4 = 0

        # c-terminus
        if length:
            if not r.l6:
                r.l6 = 0

        if not r.phi:
            r.phi = math.radians(999)
        if not r.psi:
            r.psi = math.radians(999)
        if torsion:
            if not r.ome:
                r.ome = math.radians(999)

        r.props = {
            'name': r.res_name,
            'id': r.fragment_id,
            'chain': r.chain_id,
            'phi': math.degrees(r.phi),
            'psi': math.degrees(r.psi),
                    }
        if angle:
            try:
                r.props['a1'] = math.degrees(r.a1)
                r.props['a2'] = math.degrees(r.a2)
                r.props['a3'] = math.degrees(r.a3)
                r.props['a4'] = math.degrees(r.a4)
                r.props['a5'] = math.degrees(r.a5)
                r.props['a6'] = math.degrees(r.a6)
                r.props['a7'] = math.degrees(r.a7)
            # Not a number, so probably didn't get set.
            except (TypeError, AttributeError):
                pass
        if length:
            r.props['l1'] = r.l1
            r.props['l2'] = r.l2
            r.props['l3'] = r.l3
            r.props['l4'] = r.l4
            r.props['l5'] = r.l5
            r.props['l6'] = r.l6
            r.props['l7'] = r.l7
        if torsion:
            r.props['ome'] = math.degrees(r.ome)
 
def optparse_setup():
    import optparse

    usage = """usage: %prog [options] [<geometry type> <PDB>]"""

    parser = optparse.OptionParser(version='%prog ' + version)
    parser.set_usage(usage)
    parser.set_defaults(
        verbose=False,
        compare_eh=False,
        compare_cdecg=False,
        )
    parser.add_option( \
        '-v', '--verbose', \
            action='store_true', \
            dest='verbose', \
            help='Use verbose output; defaults to %default')
    parser.add_option( \
        '-e', '--engh-and-huber', \
            action='store_true', \
            dest='compare_eh', \
            help='Compare geometry to Engh & Huber; defaults to %default')
    parser.add_option( \
        '-c', '--conformation-dependent-geometry', \
            action='store_true', \
            dest='compare_cdecg', \
            help='Compare geometry to conformation-dependent expected covalent geometry; defaults to %default')
    parser.add_option( \
        '-s', '--pdb-structure', \
            action='append', \
            dest='compare_pdbs', \
            help='Compare geometry to structure file PDB',
            metavar='PDB',
        )

    return parser

if __name__ == '__main__':
	sys.exit(main(sys.argv))
