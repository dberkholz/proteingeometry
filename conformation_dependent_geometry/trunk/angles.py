#!/usr/bin/env python

"""
angles

Purpose: returns the tau3 angle from our Protein Geometry Database,
 given a residue type and phi and psi angles.

Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
 P. Andrew Karplus laboratory, Oregon State University

Design:
 Implement as 2D dictionary matrices of bins. Bins will be classes holding
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

 If we don't have known database values, we will return the template Engh &
 Huber values.
"""

import sys, getopt

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

verbose = False

class Usage(Exception):
    """This class gives a usage message"""

    def __init__(self):
        saveout = sys.stdout
        sys.stdout = sys.stderr
        print "angles [-v|--verbose] <residue> <phi> <psi>"
        print
        print "Unknown residue names will use the 'all' residue class."
        print "The 'default' residue class is Engh & Huber values."
        print "Available residue classes:",
        for db in databases:
            print db,
        print
        sys.stdout = saveout


class bin(object):
    """This class holds all the info about a bin"""

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
        # separate on space
        words = line.split()

        # Create sorted list of available geometry attributes
        if line[0:8] == 'PhiStart':
            bin.__slots__= [slot for slot in words]
            continue

        # grab variables we need to set up the class
        phi = int(words[0])
        psi = int(words[2])

        # instantiate the class and give it a short name
        dict[(phi,psi)] = bin()
        db = dict[(phi,psi)]

        # assign values to bin
        db.PhiStart = phi
        db.PhiStop = int(words[1])
        db.PsiStart = psi
        db.PsiStop = int(words[3])
        db.Observations = int(words[4])
        db.PhiAvg = float(words[5])
        db.PhiDev = float(words[6])
        db.PsiAvg = float(words[7])
        db.PsiDev = float(words[8])
        db.L1Avg = float(words[9])
        db.L1Dev = float(words[10])
        db.L2Avg = float(words[11])
        db.L2Dev = float(words[12])
        db.L3Avg = float(words[13])
        db.L3Dev = float(words[14])
        db.L4Avg = float(words[15])
        db.L4Dev = float(words[16])
        db.L5Avg = float(words[17])
        db.L5Dev = float(words[18])
        db.a1Avg = float(words[19])
        db.a1Dev = float(words[20])
        db.a2Avg = float(words[21])
        db.a2Dev = float(words[22])
        db.a3Avg = float(words[23])
        db.a3Dev = float(words[24])
        db.a4Avg = float(words[25])
        db.a4Dev = float(words[26])
        db.a5Avg = float(words[27])
        db.a5Dev = float(words[28])
        db.a6Avg = float(words[29])
        db.a6Dev = float(words[30])
        db.a7Avg = float(words[31])
        db.a7Dev = float(words[32])
        db.OmeAvg = float(words[33])
        db.OmeDev = float(words[34])
        db.Color = float(words[35])
        db.X = float(words[36])
        db.Y = float(words[37])
        db.XLen = float(words[38])
        db.YLen = float(words[39])
        db.XBin = float(words[40])
        db.YBin = float(words[41])

    file.close()
    return dict

def get_binsize(dict):
    """Find bin size for a given database by checking a single bin in it"""
    for i in dict:
        j=dict[i]
        phi_binsize = j.PhiStop - j.PhiStart
        psi_binsize = j.PsiStop - j.PsiStart
        break
    return phi_binsize, psi_binsize

def vprint(*args):
    if verbose:
        for arg in args:
            print >> sys.stderr, arg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            optlist, args = getopt.getopt(argv[1:], 'v', ['verbose'])
        except getopt.GetoptError:
            raise Usage
        for o, a in optlist:
            if o == "-v" or o == "--verbose":
                global verbose
                verbose = True

        if len(args) != 3:
            raise Usage
        else:
            residue = args[0].lower()
            # int() can't convert from string and from float at the same time
            phi = int(float(args[1]))
            psi = int(float(args[2]))

        # Build up databases
        dblist = {}
        for database in databases:
            vprint("Creating database " + database \
                  + " with file " + databases[database])
            dblist[database] = create_database(databases[database])

        # Decide which database to use
        for database in databases:
            if residue == database:
                dbname=database
            else:
                dbname='all'

        try:
#            vprint("Database list = ", dblist)
            vprint("Database name = " + dbname)
#            vprint("Database = ", dblist[dbname])

            # Print field names
            if verbose:
                for i in dblist[dbname]:
                    j=dblist[dbname][i]
                    print '\t'.join([slot for slot in j.__slots__])
                    break

            phi_binsize, psi_binsize = get_binsize(dblist[dbname])
            print dblist[dbname][(phi-phi%phi_binsize, psi-psi%psi_binsize)]
        except KeyError:
            vprint("Defaulting to library value")

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
    except Usage:
        return 2

if __name__ == '__main__':
    # run as a program
    sys.exit(main())
