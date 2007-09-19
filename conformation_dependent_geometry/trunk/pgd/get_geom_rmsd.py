#!/usr/bin/env python

import math, sys, mmLib.Structure, mmLib.FileIO, mmLib.AtomMath
import pgd.angles as a

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
        self.dbdict = a.create_all_databases(a.databases)

    def get(self, phi, psi, attr):
        # Get based on residue type, phi and psi
        #print phi, psi
        #print type(attr)

        # angle/length-naming glue
        if attr[0] != 'a':
            pgdattr = attr.capitalize()
        else:
            pgdattr = attr
        pgdattr += 'Avg'

        if phi > 180 or psi > 180:
            residue = 'default'
            phi = 0
            psi = 0
        # Note next_res_name
        elif self.next_res_name == 'PRO':
            residue = 'glycine'
        elif self.res_name == 'GLY':
            residue = 'glycine'
        elif self.res_name == 'PRO':
            residue = 'proline'
        elif self.res_name == 'ILE' \
                or self.res_name == 'VAL':
            residue = 'ileval'
        else:
            residue = 'all'
        fields, geometry = a.get_geometry(self.dbdict, residue, phi, psi)
        #print geometry.__dict__
        value = getattr(geometry, pgdattr)
        # Our defaults don't include omega for some reason,
        # instead they return -1
        if attr == 'ome' and value == -1:
            return 180
        else:
            return value


def main(argv):
    if len(argv[1:]) != 3:
        print "Requires a length/angle and 2 PDB files as an argument"
        sys.exit(1)

    meas = argv[1]
    pdb1 = argv[2]
    pdb2 = argv[3]
    struct1 = mmLib.FileIO.LoadStructure(file=pdb1)
    struct2 = mmLib.FileIO.LoadStructure(file=pdb2)
    get_geometry(struct1)
    get_geometry(struct2)

    eh = geom(meas)
    pgd = protein_geometry_database(meas)
    #print eh.__dict__
    if meas == 'ome':
        #print "Measurement is ome"
        eh.set(meas, 180)
    msd = 0
    eh.msd = 0
    pgd.msd = 0
    N = 0
    for r in struct1.iter_amino_acids():
        r_atom = r.get_atom('C')
        r2_atom = struct2.get_equivalent_atom(r_atom)
        r2 = r2_atom.get_fragment()
        if not r2:
            continue
        if r2_atom.occupancy < 0.5 or r_atom.occupancy < 0.5:
            continue

        # We use this to look up conformation-dependent geometry
        pgd.res_name = r2.res_name
        # Look for residue i+1 for a proline, so we know if we're XPro
        pgd.next_res_numstr = str(int(r2.fragment_id) + 1)
        r2_dict = r2.get_chain().fragment_dict
        try:
            pgd.next_res_name = r2_dict[pgd.next_res_numstr]
        except:
            pgd.next_res_name = 'END'

        if meas == 'a3':
            #print "Measurement is a3"
            if r2.res_name == 'GLY':
                #print r2
                eh.set(meas, 112.5)
            elif r2.res_name == 'PRO':
                #print r2
                eh.set(meas, 111.8)
            else:
                eh.set(meas, 111.2)
        elif meas == 'phi' or meas == 'psi' or meas == 'ome':
            #print "phi/psi/ome"
            if r.props[meas] > 180 or r2.props[meas] > 180:
                continue

        dev = r.props[meas] - r2.props[meas]
        msd += dev**2
        eh.dev = r2.props[meas] - eh.get(meas)
        eh.msd += eh.dev**2
        pgd_meas = pgd.get(r.props['phi'], r.props['psi'], meas)
        pgd.dev = r2.props[meas] - pgd_meas
        pgd.msd += pgd.dev**2
        N += 1
        #print '%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d' % (r.props[meas], r2.props[meas], eh.get(meas), pgd_meas, dev, msd, eh.dev, eh.msd, pgd.dev, pgd.msd, N)
    rmsd = math.sqrt ( msd / N )
    eh.rmsd = math.sqrt ( eh.msd / N )
    pgd.rmsd = math.sqrt ( pgd.msd / N )
    print '%s for %s vs %s = %.2f' % (meas, pdb1, pdb2, rmsd)
    print '%s for %s vs E&H = %.2f' % (meas, pdb2, eh.rmsd)
    print '%s for %s vs PGD = %.2f' % (meas, pdb2, pgd.rmsd)

def get_geometry(struct):
    for r in struct.iter_amino_acids():
        next_res = r.get_offset_residue(1)
        prev_res = r.get_offset_residue(-1)

        r.a3, r.a2, r.a4, r.a5, a6, a1 = \
              r.calc_mainchain_bond_angle()
        r.l2, r.l4, r.l5, r.l3, r.l6 = \
              r.calc_mainchain_bond_length()
        r.l7 = 0
        if prev_res:
            r.l1 = prev_res.l6
        else:
            # n-terminus
            r.l1 = 0
            r.a1 = 0
            r.a6 = 0
            r.a7 = 0
        if next_res:
            next_res.a1 = a1
            next_res.a6 = a6
            aN = r.get_atom('N')
            aCA = r.get_atom('CA')
            aO = r.get_atom('O')
            aC = r.get_atom('C')
            naN = next_res.get_atom('N')

            next_res.l2 = mmLib.AtomMath.calc_distance(aN, aCA)
            r.l7 = next_res.l2
            next_res.a7 = mmLib.AtomMath.calc_angle(aO, aC, naN)
        else:
            # c-terminus
            r.l7 = 0
        r.phi = r.calc_torsion_phi()
        r.psi = r.calc_torsion_psi()
        r.ome = r.calc_torsion_omega()

        if not r.l3:
            # glycine
            r.l3 = 0
            r.a2 = 0
            r.a4 = 0
        if not r.l6:
            # c-terminus
            r.l6 = 0
        if not r.phi:
            r.phi = math.radians(999)
        if not r.psi:
            r.psi = math.radians(999)
        if not r.ome:
            r.ome = math.radians(999)

        r.props = {
            'name': r.res_name,
            'id': r.fragment_id,
            'phi': math.degrees(r.phi),
            'psi': math.degrees(r.psi),
            'ome': math.degrees(r.ome),
            'l1': r.l1,
            'l2': r.l2,
            'l3': r.l3,
            'l4': r.l4,
            'l5': r.l5,
            'l6': r.l6,
            'l7': r.l7,
            'a1': math.degrees(r.a1),
            'a2': math.degrees(r.a2),
            'a3': math.degrees(r.a3),
            'a4': math.degrees(r.a4),
            'a5': math.degrees(r.a5),
            'a6': math.degrees(r.a6),
            'a7': math.degrees(r.a7),
                    }

if __name__ == '__main__':
	sys.exit(main(sys.argv))
