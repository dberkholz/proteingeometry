#!/usr/bin/env python

# Design: Read Dunbrack file, convert into angles dictionary, then dump angles
# dictionary to file using angles.py dump feature

import bz2, sys
# Do it this way so __init__.py runs and initializes angles.optlist. Otherwise
# we get tracebacks.
from conformation_dependent_geometry.angles import bin, databases, create_all_databases

results_file = 'data/KernRegr_CDL.txt'

class_map = {
    ('NonPGIV_nonxpro', 'B'): 'other',
    ('NonPGIV_nonxpro', 'I'): 'default',
    ('Gly_nonxpro', 'B'): 'glycine',
    ('Gly_nonxpro', 'I'): 'default-gly',
    ('Pro_nonxpro', 'B'): 'proline',
    ('Pro_nonxpro', 'I'): 'default-pro',
    ('IleVal_nonxpro', 'B'): 'ileval',
    ('IleVal_nonxpro', 'I'): 'default-ileval',
    ('NonPGIV_xpro', 'B'): 'preproline',
    ('NonPGIV_xpro', 'I'): 'default-preproline',
    ('Gly_xpro', 'B'): 'glycine-preproline',
    ('Gly_xpro', 'I'): 'default-gly-preproline',
    ('Pro_xpro', 'B'): 'proline-preproline',
    ('Pro_xpro', 'I'): 'default-pro-preproline',
    ('IleVal_xpro', 'B'): 'ileval-preproline',
    ('IleVal_xpro', 'I'): 'default-ileval-preproline',
    }

bin_map = {
    'ResTypeGroup': '',
    'Phi': 'Phi',
    'Psi': 'Psi',
    'S':   '',
    'Num': 'Observations',
    }

geometry_map = {
    'CNA': 'a1',
    'NAB': 'a2',
    'NAC': 'a3',
    'BAC': 'a4',
    'ACO': 'a5',
    'ACN': 'a6',
    'OCN': 'a7',
    'CN':  'L1',
    'NA':  'L2',
    'AB':  'L3',
    'AC':  'L4',
    'CO':  'L5',
    }

type_map = {
    'm': 'Avg',
    's': 'Dev',
    }


class dunbrack_bin(bin):
    def __init__(self, words):
        #print self.var_order
        super(dunbrack_bin, self).__init__(words)


def read_dunbrack(infile):
    dunbrack_class_dict = {}

    for line in open(infile):
        if line[0] == '#':
            continue
        # Descriptive information
        elif line[0] == '@':
            line = line[1:]
        words = line.split()
        if not words:
            continue
        elif words[0] == 'Format':
            format = words[1]
            continue
        elif words[0] == 'Version':
            version = words[1]
            continue

        # Create sorted list of available geometry attributes
        if words[0] == 'ResTypeGroup':
            dunbrack_bin.var_order = words
            continue

        # grab variables we need to set up the class
        phi = int(words[1])
        psi = int(words[2])
        dunbrack_class = words[0]

        try:
            dunbrack_class_dict[dunbrack_class]
        except KeyError:
            dunbrack_class_dict[dunbrack_class] = {}

        if (phi,psi) in dunbrack_class_dict[dunbrack_class]:
            print 'Overwriting', phi,psi
        dunbrack_class_dict[dunbrack_class][(phi,psi)] = dunbrack_bin(words)

    #print dunbrack_class_dict.keys()
    #for key in dunbrack_class_dict:
    #    print dunbrack_class_dict[key][(-40,-60)]
    #    print
    return dunbrack_class_dict


def translate_dunbrack_bin(dunbrack_bin):
    database_bin = bin(words='')

    # We also need to set up PhiStart, PhiStop, PsiStart, PsiStop,
    # Observations.  The file is missing omega, chi, zeta, hbond.
    for bin_name in bin_map:
        dunbrack_field = bin_name
        if not hasattr(dunbrack_bin, dunbrack_field):
            continue
        pgd_field = bin_map[bin_name]
        if pgd_field is 'Phi' or pgd_field is 'Psi':
            pgd_start_field = pgd_field + 'Start'
            pgd_stop_field = pgd_field + 'Stop'
            # Since we use floor rounding, and Dunbrack's numbers are supplied
            # as point values rather than ranges for bins, pretend they are
            # 10-degree bins with range (value-5, value+5)
            setattr(database_bin, pgd_start_field, getattr(dunbrack_bin, dunbrack_field) - 5)
            setattr(database_bin, pgd_stop_field, getattr(dunbrack_bin, dunbrack_field) + 5)
        else:
            setattr(database_bin, pgd_field, getattr(dunbrack_bin, dunbrack_field))

    for geometry in geometry_map:
        for geom_type in type_map:
            dunbrack_field = geom_type + geometry
            pgd_field = geometry_map[geometry] + type_map[geom_type]
            if hasattr(dunbrack_bin, dunbrack_field):
                setattr(database_bin, pgd_field, getattr(dunbrack_bin, dunbrack_field))

    # Here we access the data structure and store our data in it
    #res_class = words[0]
    #database_dict[class_map[res_class]]

    return database_bin


def convert_dunbrack_database():
    # The only reason we run this is to set up the angles.bin attributes list
    create_all_databases(databases)

    dunbrack_class_dict = read_dunbrack(results_file)
    unused_attribs = ['PhiAvg', 'PhiDev', 'PsiAvg', 'PsiDev', 'OmeAvg', 'OmeDev', 'ChiAvg', 'ChiDev', 'ZetaAvg', 'ZetaDev', 'HBondAvg', 'HBondDev']
    for attrib in unused_attribs:
        bin.var_order.remove(attrib)

    # Convert data and print it to file
    for dunbrack_dict in dunbrack_class_dict:
        # We only want to print backbone-independent data once
        found_independent = False

        filename_base = 'data/1.0-dunbrack-'
        independent = bz2.BZ2File(filename_base + class_map[(dunbrack_dict, 'I')] + '.bz2', 'w')
        dependent   = bz2.BZ2File(filename_base + class_map[(dunbrack_dict, 'B')] + '.bz2', 'w')

        # Call angles.get_fields(pgd_dict) to print field names

        for dunbrack_bin in dunbrack_class_dict[dunbrack_dict]:
            pgd_bin = translate_dunbrack_bin(dunbrack_class_dict[dunbrack_dict][dunbrack_bin])
            if dunbrack_class_dict[dunbrack_dict][dunbrack_bin].S == 'I':
                if found_independent:
                    continue
                print >> independent, pgd_bin
                found_independent = True
            elif dunbrack_class_dict[dunbrack_dict][dunbrack_bin].S == 'B':
                print >> dependent, pgd_bin

        independent.close()
        dependent.close()


def main():
    convert_dunbrack_database()


if __name__ == '__main__':
    # run as a program
    sys.exit(main())
