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
 into the dictionary matrix. Dictionary keys will be tuples of phi,psi pairs.

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

databases = {
'all' : '1.0-graphdata-all-but-gly-pro-xpro.txt',
'glycine' : '1.0-graphdata-gly.txt',
'proline': '1.0-graphdata-pro.txt',
'preproline': '1.0-graphdata-xpro.txt'
}

class bin(object):
    """This class holds all the info about a bin"""

    def __init__(self):
        """initialize attributes to obviously incorrect values"""
        self.phi_start = -1
        self.phi_stop = -1
        self.psi_start = -1
        self.psi_stop = -1
        self.observations = -1
        self.phi_avg = -1
        self.phi_dev = -1
        self.psi_avg = -1
        self.psi_dev = -1
        self.l1_avg = -1
        self.l1_dev = -1
        self.l2_avg = -1
        self.l2_dev = -1
        self.l3_avg = -1
        self.l3_dev = -1
        self.l4_avg = -1
        self.l4_dev = -1
        self.l5_avg = -1
        self.l5_dev = -1
        self.a1_avg = -1
        self.a1_dev = -1
        self.a2_avg = -1
        self.a2_dev = -1
        self.a3_avg = -1
        self.a3_dev = -1
        self.a4_avg = -1
        self.a4_dev = -1
        self.a5_avg = -1
        self.a5_dev = -1
        self.a6_avg = -1
        self.a6_dev = -1
        self.a7_avg = -1
        self.a7_dev = -1
        self.omega_avg = -1
        self.omega_dev = -1

    def __str__(self):
        """Set up so we can just print class instances to get the output we
        want"""

        # Can we generalize this to just print all class attributes?
        print \
              self.phi_start, \
              self.phi_stop, \
              self.psi_start, \
              self.psi_stop, \
              self.observations, \
              self.phi_avg, \
              self.phi_dev, \
              self.psi_avg, \
              self.psi_dev, \
              self.l1_avg, \
              self.l1_dev, \
              self.l2_avg, \
              self.l2_dev, \
              self.l3_avg, \
              self.l3_dev, \
              self.l4_avg, \
              self.l4_dev, \
              self.l5_avg, \
              self.l5_dev, \
              self.a1_avg, \
              self.a1_dev, \
              self.a2_avg, \
              self.a2_dev, \
              self.a3_avg, \
              self.a3_dev, \
              self.a4_avg, \
              self.a4_dev, \
              self.a5_avg, \
              self.a5_dev, \
              self.a6_avg, \
              self.a6_dev, \
              self.a7_avg, \
              self.a7_dev, \
              self.omega_avg, \
              self.omega_dev

class default_bin(bin):
    def __init__(self):
        """initialize attributes to Engh & Huber default values"""
        self.phi_start = -1
        self.phi_stop = -1
        self.psi_start = -1
        self.psi_stop = -1
        self.observations = -1
        self.phi_avg = -1
        self.phi_dev = -1
        self.psi_avg = -1
        self.psi_dev = -1
        self.l1_avg = -1
        self.l1_dev = -1
        self.l2_avg = -1
        self.l2_dev = -1
        self.l3_avg = -1
        self.l3_dev = -1
        self.l4_avg = -1
        self.l4_dev = -1
        self.l5_avg = -1
        self.l5_dev = -1
        self.a1_avg = -1
        self.a1_dev = -1
        self.a2_avg = -1
        self.a2_dev = -1
        self.a3_avg = -1
        self.a3_dev = -1
        self.a4_avg = -1
        self.a4_dev = -1
        self.a5_avg = -1
        self.a5_dev = -1
        self.a6_avg = -1
        self.a6_dev = -1
        self.a7_avg = -1
        self.a7_dev = -1
        self.omega_avg = -1
        self.omega_dev = -1

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

    file.close()
    return dict

def usage():
    """Print usage message and quit"""
    print "angles <residue> <phi> <psi>"
    sys.exit(1)

if __name__ == '__main__':
    # run as a program
    if len(sys.argv) != 4:
        usage()
    else:
        residue = sys.argv[1]
        phi = sys.argv[2]
        psi = sys.argv[3]

        # round phi/psi down to nearest <binsize>
        phi = int(phi)/10*10
        psi = int(psi)/10*10

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
            print dblist[dbname][(phi,psi)]
        except:
            try:
                default = default_bin()
                print default
            except TypeError:
                pass
