#!/usr/bin/env python

"""
angles

Purpose: returns the tau3 angle from our Protein Geometry Database,
 given a residue type and phi and psi angles.

Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
 P. Andrew Karplus laboratory, Oregon State University

Design:
 Implement as a 2D dictionary matrix of bins. Bins will be classes holding
 all info about that bin. You can look up a bin by passing its phi,psi pair
 into the dictionary matrix. Dictionary keys will be tuples of phi,psi
 pairs.

 Databases will be read into the dictionary matrices from files. The best
 format will be one file per residue class (proline, glycine, preproline,
 all others). This gives clean conceptual separation, where residue classes
 are equivalent to files. It also ensures that each line contains a single
 bin's geometry from a single residue class to ease parsing. Each residue
 class will have its own dictionary matrix.

 Angles passed in as arguments will be rounded down to the nearest bin start
 and then passed to the matrix to return a bin class.

 We can use the has_key() dictionary method to check whether we have known
 database values, and if not, return the template of Engh & Huber values.
"""

import sys

# Database containing fallback, standard values (likely Engh & Huber)
default = 'default'

# Map of all residue classes and the file names containing their libraries
# The key names are the valid residue classes to pass in as arguments
databases = {
'all' : '1.0-graphdata-all-but-gly-pro-xpro.txt',
'glycine' : '1.0-graphdata-gly.txt',
'proline': '1.0-graphdata-pro.txt',
'preproline': '1.0-graphdata-xpro.txt',
'default': 'engh-huber.txt'
}

class bin(object):
    """This class holds all the info about a bin"""

    # Define legal attributes
    __slots__ = [
        'phi_start',
        'phi_stop',
        'psi_start',
        'psi_stop',
        'observations',
        'phi_avg',
        'phi_dev',
        'psi_avg',
        'psi_dev',
        'l1_avg',
        'l1_dev',
        'l2_avg',
        'l2_dev',
        'l3_avg',
        'l3_dev',
        'l4_avg',
        'l4_dev',
        'l5_avg',
        'l5_dev',
        'a1_avg',
        'a1_dev',
        'a2_avg',
        'a2_dev',
        'a3_avg',
        'a3_dev',
        'a4_avg',
        'a4_dev',
        'a5_avg',
        'a5_dev',
        'a6_avg',
        'a6_dev',
        'a7_avg',
        'a7_dev',
        'omega_avg',
        'omega_dev'
        ]

    def __str__(self):
        """Print all class attributes"""

        list = [str(getattr(self, attr)) for attr in self.__slots__]
        return '\t'.join(list)

def create_database(filename):
    """Create a dictionary matrix holding all of the bins

    Returns the dictionary matrix"""
    file = open(filename, 'r')
    dict = {}
    for line in file:
        # Don't read the first line
        if line[0:8] == 'PhiStart':
            continue

        # separate on space
        words = line.split()

        # grab variables we need to set up the class
        phi = int(words[0])
        psi = int(words[2])

        # instantiate the class and give it a short name
        dict[(phi,psi)] = bin()
        db = dict[(phi,psi)]

        # assign values to bin
        db.phi_start = phi
        db.phi_stop = int(words[1])
        db.psi_start = psi
        db.psi_stop = int(words[3])
        db.observations = int(words[4])
        db.phi_avg = float(words[5])
        db.phi_dev = float(words[6])
        db.psi_avg = float(words[7])
        db.psi_dev = float(words[8])
        db.l1_avg = float(words[9])
        db.l1_dev = float(words[10])
        db.l2_avg = float(words[11])
        db.l2_dev = float(words[12])
        db.l3_avg = float(words[13])
        db.l3_dev = float(words[14])
        db.l4_avg = float(words[15])
        db.l4_dev = float(words[16])
        db.l5_avg = float(words[17])
        db.l5_dev = float(words[18])
        db.a1_avg = float(words[19])
        db.a1_dev = float(words[20])
        db.a2_avg = float(words[21])
        db.a2_dev = float(words[22])
        db.a3_avg = float(words[23])
        db.a3_dev = float(words[24])
        db.a4_avg = float(words[25])
        db.a4_dev = float(words[26])
        db.a5_avg = float(words[27])
        db.a5_dev = float(words[28])
        db.a6_avg = float(words[29])
        db.a6_dev = float(words[30])
        db.a7_avg = float(words[31])
        db.a7_dev = float(words[32])
        db.omega_avg = float(words[33])
        db.omega_dev = float(words[34])

    file.close()
    return dict

def get_binsize(dict):
    """Find bin size for a given database by checking a single bin in it"""
    for i in dict:
        j=dict[i]
        phi_binsize = j.phi_stop - j.phi_start
        psi_binsize = j.psi_stop - j.psi_start
        break
    return phi_binsize, psi_binsize

def usage():
    """Print usage message and quit"""
    print "angles <residue> <phi> <psi>"
    print
    print "Unknown residue names will use the 'all' residue class."
    print "The 'default' residue class is Engh & Huber values."
    print "Available residue classes:",
    for db in databases:
        print db,
    sys.exit(1)

if __name__ == '__main__':
    # run as a program
    if len(sys.argv) != 4:
        usage()
    else:
        residue = sys.argv[1].lower()
        # int() can't convert from string and from float at the same time
        phi = int(float(sys.argv[2]))
        psi = int(float(sys.argv[3]))

        # Build up databases
        dblist = {}
        for database in databases:
            print "Creating database " + database \
                  + " with file " + databases[database]
            dblist[database] = create_database(databases[database])

        # Decide which database to use
        for database in databases:
            if residue == database:
                dbname=database
            else:
                dbname='all'

        try:
            print
#            print "dblist = ", dblist
            print "dbname = " + dbname
#            print "phi, psi = ", phi, psi
#            print dblist[dbname]
            phi_binsize, psi_binsize = get_binsize(dblist[dbname])
            print dblist[dbname][(phi-phi%phi_binsize, psi-psi%psi_binsize)]
        except KeyError:
            print "Defaulting to library value"

            # Special code to deal with huge binsizes
            # e.g., 360x360 as in the default library
            phi_binsize, psi_binsize = get_binsize(dblist[default])
            phi_div, phi_mod = divmod(abs(phi), phi_binsize)
            psi_div, psi_mod = divmod(abs(psi), psi_binsize)
            if phi_div == 0:
                phi_r = 180 - phi_binsize
            else:
                phi_r = phi-phi_mod
            if psi_div == 0:
                psi_r = 180 - psi_binsize
            else:
                psi_r = psi-psi_mod

            print dblist[default][(phi_r, psi_r)]
