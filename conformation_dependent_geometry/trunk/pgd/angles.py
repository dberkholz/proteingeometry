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

import os, sys, optparse

# Program version
version = '0.1'

# Database containing fallback, standard values (likely Engh & Huber)
default = 'default'

# Map of all residue classes and the file names containing their libraries
# The key names are the valid residue classes to pass in as arguments
moduledir = os.path.dirname(__file__)
databases = {
'all' : moduledir + '/data/1.0-graphdata-all-but-gly-pro-xpro.txt',
'glycine' : moduledir + '/data/1.0-graphdata-gly.txt',
'proline': moduledir + '/data/1.0-graphdata-pro.txt',
'preproline': moduledir + '/data/1.0-graphdata-xpro.txt',
'default': moduledir + '/data/engh-huber.txt'
}

# Option handling
dblist = ' '.join(databases)
usage = """usage: %prog [options] <residue> <phi> <psi>

Unknown residue names will use the 'all' residue class.
The 'default' residue class is Engh & Huber values.
Available residue classes: """ + dblist

parser = optparse.OptionParser(version='%prog ' + version)
parser.disable_interspersed_args()
parser.set_usage(usage)
parser.set_defaults(verbose=False)
parser.add_option( \
    '-v', '--verbose', \
        action='store_true', \
        dest='verbose', \
        help='Use verbose output; defaults to %default')
optlist, args = parser.parse_args()

class bin(object):
    """This class holds all the info about a bin"""

    def __str__(self):
        """Print all class attributes"""

        return '\t'.join(str(getattr(self, attr)) for attr in self.var_order)

def create_database(filename):
    """Create a dictionary matrix holding all of the bins

    Returns the dictionary matrix"""
    dbdict = {}
    for line in open(filename, 'r'):
        # separate on space
        words = line.split()

        # Create sorted list of available geometry attributes
        if line.startswith('PhiStart'):
            bin.var_order = words
            continue

        # grab variables we need to set up the class
        phi = int(words[0])
        psi = int(words[2])

        # instantiate the class and give it a short name
        dbdict[(phi,psi)] = bin()
        db = dbdict[(phi,psi)]

        # assign values to bin
        for i, slot in enumerate(db.var_order):
            if '.' in words[i]:
                words[i] = float(words[i])
            else:
                words[i] = int(words[i])
            setattr(db, slot, words[i])

    return dbdict

def get_binsize(dbdict):
    """Find bin size for a given database by checking a single bin in it"""
    for v in dbdict.itervalues():
        phi_binsize = v.PhiStop - v.PhiStart
        psi_binsize = v.PsiStop - v.PsiStart
        return phi_binsize, psi_binsize

def get_default_binsize(phi, psi, phi_binsize, psi_binsize):
    """Find the default bin size, using specialized code for large bins"""

    phi_div, phi_mod = divmod(abs(phi), phi_binsize)
    psi_div, psi_mod = divmod(abs(psi), psi_binsize)
    if not phi_div:
        phi_r = 180 - phi_binsize
    else:
        phi_r = phi-phi_mod
    if not psi_div:
        psi_r = 180 - psi_binsize
    else:
        psi_r = psi-psi_mod
    return phi_r, psi_r

def get_fields(database):
    """Return field names"""
    if optlist.verbose:
        database.keys().sort()
        for v in database.itervalues():
            fields = '\t'.join(v.var_order)
            return fields

def get_geometry(dblist, residue, phi, psi):
    """Get the geometry info for a specific residue/phi/psi setting

    Return the field names and the geometry."""

    # Decide which database to use
    for database in databases:
        if residue == database:
            dbname=database
        else:
            dbname='all'

    fields = get_fields(dblist[dbname])

    try:
        #vprint("Database list = ", dblist)
        vprint("Database name = " + dbname)
        #vprint("Database = ", dblist[dbname])

        phi_binsize, psi_binsize = get_binsize(dblist[dbname])
        return fields, dblist[dbname][(phi-phi%phi_binsize, psi-psi%psi_binsize)]
    except KeyError:
        vprint("Defaulting to library value")
        dbname=default

        phi_binsize, psi_binsize = get_binsize(dblist[dbname])
        phi_r, psi_r = get_default_binsize(phi, psi, phi_binsize, psi_binsize)
        return fields, dblist[dbname][(phi_r, psi_r)]

def vprint(*args):
    """Verbose print; only print if verbosity is enabled."""

    if optlist.verbose:
        for arg in args:
            print >> sys.stderr, arg

def main(argv):
    if len(args) != 3:
        parser.error('incorrect number of arguments')
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

    fields, geometry = get_geometry(dblist, residue, phi, psi)
    if fields:
        print fields
    print geometry

if __name__ == '__main__':
    # run as a program
    sys.exit(main(sys.argv))
