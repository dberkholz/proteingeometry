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

"""
angles

Purpose: returns geometry values from our Protein Geometry Database,
 given a residue, the next residue, and phi and psi angles.

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

import math, os, sys, bz2

# Program version
version = '0.3'

# Database containing fallback, standard values (likely Engh & Huber)
default = 'default'

# Minimum number of observations required to return non-E&H geometry
observation_min = 3

# Map of all residue classes and the file names containing their libraries
# The key names are the valid residue classes to pass in as arguments
moduledir = os.path.dirname(__file__)

# General design of using Dunbrack data: A dictionary of bins should first
# look to Dunbrack, then automatically and transparently fall back to
# PGD. This will require knowledge of backup dictionaries at the dbdict level,
# and it needs to be smart enough to do this: Dunbrack backbone-dependent, PGD
# backbone-dependent, Dunbrack global average, PGD global average. So we will
# need a new class that holds two dictionaries in parallel per residue type
# and does the magic when a bin is requested.
backbone_independent_databases = {}
databases = {
'independent': {},
'dependent': {},
}

databases['independent']['dunbrack'] = {
    'default' : moduledir + '/data/1.0-dunbrack-default.txt.bz2',
    'glycine' : moduledir + '/data/1.0-dunbrack-default-glycine.txt.bz2',
    'ileval': moduledir + '/data/1.0-dunbrack-default-ileval.txt.bz2',
    'proline': moduledir + '/data/1.0-dunbrack-default-proline.txt.bz2',
    'preproline': moduledir + '/data/1.0-dunbrack-default-preproline.txt.bz2',
    ('glycine','preproline') : moduledir + '/data/1.0-dunbrack-default-glycine-preproline.txt.bz2',
    ('ileval','preproline'): moduledir + '/data/1.0-dunbrack-default-ileval-preproline.txt.bz2',
    ('proline','preproline'): moduledir + '/data/1.0-dunbrack-default-proline-preproline.txt.bz2',
}

databases['independent']['eh'] = {
    'default': moduledir + '/data/engh-huber.txt.bz2',
    'glycine': moduledir + '/data/engh-huber-gly.txt.bz2',
    'proline': moduledir + '/data/engh-huber-pro.txt.bz2',
}

databases['dependent']['dunbrack'] = {
    'other' : moduledir + '/data/1.0-dunbrack-other.txt.bz2',
    'glycine' : moduledir + '/data/1.0-dunbrack-glycine.txt.bz2',
    'ileval': moduledir + '/data/1.0-dunbrack-ileval.txt.bz2',
    'proline': moduledir + '/data/1.0-dunbrack-proline.txt.bz2',
    'preproline': moduledir + '/data/1.0-dunbrack-preproline.txt.bz2',
    ('glycine','preproline') : moduledir + '/data/1.0-dunbrack-glycine-preproline.txt.bz2',
    ('ileval','preproline'): moduledir + '/data/1.0-dunbrack-ileval-preproline.txt.bz2',
    ('proline','preproline'): moduledir + '/data/1.0-dunbrack-proline-preproline.txt.bz2',
}

databases['dependent']['pgd'] = {
    'other' : moduledir + '/data/1.0-graphdata-all-but-gly-pro-xpro-ile-val.txt.bz2',
    'glycine' : moduledir + '/data/1.0-graphdata-gly.txt.bz2',
    'ileval': moduledir + '/data/1.0-graphdata-ile-val.txt.bz2',
    'proline': moduledir + '/data/1.0-graphdata-pro.txt.bz2',
    'preproline': moduledir + '/data/1.0-graphdata-xpro.txt.bz2',
}

# Angles are defined as tuples. The number indicates which residue of a
# 3-residue sequence a given atom is in

# These names are keys into the expected_atom_i_seqs
# Order: phi, psi, omega
dihedral_names = [ 'phi', 'psi', 'ome' ]
dihedral_atoms = (
    ((0, 'C' ), (1, 'N' ), (1, 'CA'), (1, 'C' )),
    ((1, 'N' ), (1, 'CA'), (1, 'C' ), (2, 'N' )),
    ((1, 'CA'), (1, 'C' ), (2, 'N' ), (2, 'CA')),
    )

# These names are keys into the expected_atom_i_seqs
# Order: a1-a7, defined in Karplus 1996
angle_names = [ 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7' ]
angle_atoms = (
    ((0, 'C' ), (1, 'N' ), (1, 'CA')),
    ((1, 'N' ), (1, 'CA'), (1, 'CB')),
    ((1, 'N' ), (1, 'CA'), (1, 'C' )),
    ((1, 'CB'), (1, 'CA'), (1, 'C' )),
    ((1, 'CA'), (1, 'C' ), (1, 'O' )),
    ((1, 'CA'), (1, 'C' ), (2, 'N' )),
    ((1, 'O' ), (1, 'C' ), (2, 'N' )),
    )

def get_database_attribute_average_name(geometry_name):
    """Accepts a name of an angle or length and returns the name of
    that attribute's average in the database"""
    if geometry_name[0] != 'a':
        average_name = geometry_name.capitalize()
    else:
        average_name = geometry_name
    average_name += 'Avg'

    return average_name

def get_database_attribute_deviation_name(geometry_name):
    """Accepts a name of an angle or length and returns the name of
    that attribute's deviation in the database"""
    if geometry_name[0] != 'a':
        deviation_name = geometry_name.capitalize()
    else:
        deviation_name = geometry_name
    deviation_name += 'Dev'

    return deviation_name

def get_residue_type(residue, next_residue):
    """Accepts two residue names in standard PDB format and returns the
    residue type the first residue belongs to. Used internally by
    get_geometry()."""

    # Types
    glycine = ['GLY']
    proline = ['PRO']
    ileval = ['ILE', 'VAL']
    other = [
        'ALA', 
        'ARG', 
        'ASN', 
        'ASP', 
        'CYS', 
        'GLN', 
        'GLU', 
        'LEU', 
        'LYS', 
        'MET', 
        'PHE', 
        'SER', 
        'THR', 
        'TRP', 
        'TYR', 
        'VAL'
        ]

    residue_type = []
    if residue in glycine:
        residue_type.append('glycine')
    elif residue in proline:
        residue_type.append('proline')
    elif residue in ileval:
        residue_type.append('ileval')
    elif residue in other:
        residue_type.append('other')

    if next_residue in proline:
        residue_type.append('preproline')

    vprint('type =', residue_type)

    return residue_type

class bin(object):
    """This class holds all the info about a bin. Bins cannot be instantiated
    without first filling in var_order. See create_database() for that code."""

    def __init__(self, words):
        """Fill in a bin, reading geometry info from words"""
        self.local_var_order = self.var_order

        if not len(words):
            return

        # assign values to bin
        for i, slot in enumerate(self.local_var_order):
            try:
                if '.' in words[i]:
                    words[i] = float(words[i])
                else:
                    words[i] = int(words[i])
            # Handle as a string
            except ValueError:
                pass
            setattr(self, slot, words[i])

    def __str__(self):
        """Print all class attributes"""

        ret_list = []
        for attr in self.local_var_order:
            try:
                ret_list.append( str( getattr(self, attr) ) )
            except AttributeError:
                continue
        return '\t'.join(ret_list)

class InvalidSource(BaseException):
    def __init__(self):
        BaseException.__init__()

class geometry_getter(object):
    def __init__(self):
        # All the data will get stored here
        self.sources = {}

    # Load dictionaries into object
    def load(self, dblist, source='pgd', method='dependent'):
        try:
            self.sources[method]
        except:
            self.sources[method] = {}
        self.sources[method][source] = dblist

    # We'll use dictionary lookups to access geometry info
    def __getitem__(self, in_tuple):
        # Receive phi/psi/residue info, find the right dictionary, return the
        # info.
        residue, next_residue, phi, psi = in_tuple

        vprint("residue = " + residue)

        residue_type = get_residue_type(residue, next_residue)

        # Priority order for methods. We'd rather fallback to Dunbrack's
        # global average than use PGD conformation-dependent numbers. PGD
        # numbers are only for parameters Dunbrack doesn't have yet.
        for source_name in 'dunbrack', 'pgd', 'eh':
            for method in 'dependent', 'independent':
                vprint("Source = " + source_name)
                vprint("Method = " + method)
                # Decide which residue class to use
                residue_class = None
                try:
                    source = self.sources[method][source_name]
                except:
                    # That source-method combination don't exist
                    continue
                residue_classes = source.keys()
                residue_class = [i for i in residue_classes \
                                 if tuple(residue_type) == i]
                if not residue_class:
                    if 'preproline' in residue_type:
                        residue_class = ['preproline']
                    if not residue_class:
                        residue_class = [i for i in residue_classes \
                                         if residue_type[0] == i]
                        if not residue_class:
                            residue_class = ['default']

                residue_class = ''.join(residue_class)

                res_db = source[residue_class]
                # Do the lookup
                fields = get_fields(res_db)


                try:
                    vprint("Database name = " + residue_class)
                    phi_binsize, psi_binsize = get_binsize(res_db)
                    vprint("binsizes:", phi_binsize, psi_binsize)
                    # find closest key -- NOTE: this will find a key no matter
                    # what, even if it's across the Ramachandran plot
                    dist, r_phi, r_psi = min((math.sqrt((phi-i)**2+(psi-j)**2),i,j)  for i,j in res_db.keys())
                    # If the distance is too far, we didn't find a key
                    if method == 'dependent':
                        if dist > math.sqrt((phi_binsize/2)**2+(psi_binsize/2)**2):
                            vprint('Skipping', dist, r_phi, r_psi)
                            continue                    
                    # If the key doesn't exist, you get a weird error:
                    # TypeError: unhashable type: 'dict'
                    geometry = res_db[(r_phi, r_psi)]
                    if getattr(geometry, 'Observations') < observation_min:
                        raise KeyError
                except KeyError:
                    # Try the next library
                    continue

                if geometry:
                    return fields, geometry


def create_database(filename):
    """Create a dictionary matrix holding all of the bins

    Returns the dictionary matrix"""
    dbdict = {}
    for line in bz2.BZ2File(filename, 'r'):
        # separate on space
        words = line.split()

        # Create sorted list of available geometry attributes
        if line.startswith('PhiStart'):
            bin.var_order = words
            continue

        # grab variables we need to set up the class
        phi = int(float(words[0]))
        psi = int(float(words[2]))

        # instantiate the class and fill it
        dbdict[(phi,psi)] = bin(words)

    return dbdict

def create_all_databases(db_names):
    """Build up databases"""
    dbdict = {}
    for database in db_names:
        vprint("Creating database " + str(database) \
                   + " with file " + db_names[database])
        dbdict[database] = create_database(db_names[database])
    return dbdict

def get_binsize(dbdict):
    """Find bin size for a given database by checking a single bin in it"""
    for v in dbdict.itervalues():
        phi_binsize = v.PhiStop - v.PhiStart
        psi_binsize = v.PsiStop - v.PsiStart
        return int(phi_binsize), int(psi_binsize)

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
            fields = '\t'.join(v.local_var_order)
            return fields

def get_database_name(residue, next_residue):
    vprint("residue = " + residue)

    residue_type = get_residue_type(residue, next_residue)

    # Decide which database to use
    for database in databases:
        #vprint("database = " + database)
        if residue_type == database:
            dbname = database
            break
        else:
            dbname = 'other'

    return dbname

def jackknife_geometry(dblist, residue, next_residue, phi, psi, param, value, jackknife_deviation=False):
    dbname = get_database_name(residue, next_residue)
    try:
        vprint("Database name = " + dbname)
        phi_binsize, psi_binsize = get_binsize(dblist[dbname])
        vprint("binsizes:", phi_binsize, psi_binsize)

        geometry = dblist[dbname][(phi-phi%phi_binsize, psi-psi%psi_binsize)]
        old_obs = getattr(geometry, 'Observations')
        if old_obs < observation_min:
            raise KeyError

        old_avg = getattr(
            geometry,
            get_database_attribute_average_name(param))
        old_dev = getattr(
            geometry,
            get_database_attribute_deviation_name(param))

        new_obs = old_obs - 1

        # Algorithms from Knuth
        new_avg = old_avg + (old_avg - value)/new_obs

        if jackknife_deviation:

            old_S = old_dev**2 * (old_obs - 1)
        
            S_difference = (value - old_avg) * (value - new_avg)
            new_S = old_S - S_difference

            try:
                # Lack of precision in the input numbers probably explains
                # the instability
                new_dev = math.sqrt(new_S / (new_obs - 1))
                dev_difference = new_dev - old_dev
                if abs(dev_difference) > 0.2:
                    print '%.2f' % dev_difference
            except ValueError:
                print '%.2f %.2f %.2f %.2f %.2f %.2f %.2f %d' % \
                      (new_S, old_S, S_difference, old_dev, new_avg, old_avg, value, new_obs)
                raise ValueError

            setattr(
                geometry,
                get_database_attribute_deviation_name(param),
                new_dev)

        # Write the new values to the dictionary
        setattr(
            geometry,
            get_database_attribute_average_name(param),
            new_avg)
        setattr(
            geometry,
            'Observations',
            new_obs)

    # We don't change the Engh & Huber fallbacks
    except KeyError:
        pass

def get_geometry(dblist, residue, next_residue, phi, psi):
    """Get the geometry info for a specific residue/phi/psi setting

    Return the field names and the geometry."""

    dbname = get_database_name(residue, next_residue)
    fields = get_fields(dblist[dbname])

    try:
        vprint("Database name = " + dbname)
        phi_binsize, psi_binsize = get_binsize(dblist[dbname])
        vprint("binsizes:", phi_binsize, psi_binsize)

        geometry = dblist[dbname][(phi-phi%phi_binsize, psi-psi%psi_binsize)]
        if getattr(geometry, 'Observations') < observation_min:
            raise KeyError
    except KeyError:
        library = default
        if dbname == 'proline':
            library = 'default-proline'
        elif dbname == 'glycine':
            library = 'default-glycine'
        dbname = library
        vprint("Defaulting to library value " + library)

        phi_binsize, psi_binsize = get_binsize(dblist[dbname])
        phi_r, psi_r = get_default_binsize(phi, psi, phi_binsize, psi_binsize)
        geometry = dblist[dbname][(phi_r, psi_r)]
    return fields, geometry

def iterate_over_bins(dbdict):
    """Adds default geometry info to zero-observation bins"""
    # FIXME: Should we first drop back to 'other', then 'default'?

    # Iterate over all bins from -180 to +180 using binsize chunks
    phi_def_binsize, psi_def_binsize = get_binsize(dbdict[default])
    for dbname in dbdict:
        if dbname == default:
            continue
        db = dbdict[dbname]
        phi_binsize, psi_binsize = get_binsize(db)
        for phi in xrange(-180, +181, phi_binsize):
            for psi in xrange(-180, +181, psi_binsize):
                # Check whether a bin instance exists
                # If not, make one with defaults
                if optlist.add_empty:
                    add_empty_bins(dbdict, db, phi, psi, \
                                       phi_def_binsize, psi_def_binsize, \
                                       phi_binsize, psi_binsize)

                # Need to dump after adding empty
                if optlist.dump == dbname:
                    if (phi, psi) in db:
                        print db[(phi, psi)]

    return dbdict

def add_empty_bins(dbdict, db, phi, psi, phidefbin, psidefbin, phibin, psibin):
    if not (phi, psi) in db:
        phi_def, psi_def = get_default_binsize( \
            phi, psi, phidefbin, psidefbin)
        db[(phi, psi)] = bin([])
        db[(phi, psi)] = dbdict[default][(phi_def, psi_def)]
        # Fix phi/psi; they were the values from 'default'
        # Otherwise binsize calculation is wrong
        db[(phi, psi)].PhiStart = phi
        db[(phi, psi)].PsiStart = psi
        db[(phi, psi)].PhiStop = phi + phibin
        db[(phi, psi)].PsiStop = psi + psibin

def vprint(*args):
    """Verbose print; only print if verbosity is enabled."""

    if optlist.verbose:
        print >> sys.stderr, ' '.join([str(arg) for arg in args])

def optparse_setup():
    import optparse

    # Option handling
    dblist = ' '.join(sorted(databases))
    usage = """usage: %prog [options] [<residue> <next_residue> <phi> <psi>]

 When no observations exist in the original residue group, the 'default' group
 of Engh & Huber values is used as a fallback. This is indicated by an
 observations count of -1.

Angles:
 Both angles 'phi' and 'psi' accept either integers or floats between -180 and
 +180.

The arguments <residue>, <phi> and <psi> are required unless dumping a
database.

For definitions of the angles and lengths, refer to karplus-definitions.jpg,
which installs with the documentation."""

    parser = optparse.OptionParser(version='%prog ' + version)
    parser.disable_interspersed_args()
    parser.set_usage(usage)
    parser.set_defaults(verbose=False, add_empty=False)
    parser.add_option( \
        '-v', '--verbose', \
            action='store_true', \
            dest='verbose', \
            help='Use verbose output; defaults to %default')
    parser.add_option( \
        '-e', '--add-empty', \
            action='store_true', \
            dest='add_empty', \
            help='Add default geometry to zero-observation bins; defaults to %default')
    parser.add_option( \
        '-d', '--dump-database', \
            dest='dump', \
            help='Dump a residue class CLASS to standard output', \
            metavar='CLASS')
    return parser

def main(argv):
    global optlist
    global args
    parser = optparse_setup()
    optlist, args = parser.parse_args()

    if not optlist.dump:
        if len(args) != 4:
            parser.error('incorrect number of arguments')
        else:
            residue = args[0].upper()
            next_residue = args[1].upper()
            # int() can't convert from string and from float at the same time
            phi = int(float(args[2]))
            psi = int(float(args[3]))

    geom = geometry_getter()

    for method in databases:
        for source in databases[method]:
            dblist = create_all_databases(databases[method][source])
            geom.load(dblist=dblist,
                      source=source,
                      method=method)

    if optlist.add_empty or optlist.dump:
        if optlist.dump:
            fields = get_fields(dblist[optlist.dump])
            if fields:
                print fields
        dblist = iterate_over_bins(dblist)

    if not optlist.dump:
        vprint("phi, psi:", phi, psi)
        #fields, geometry = get_geometry(dblist, residue, next_residue, phi, psi)
        fields, geometry = geom[residue, next_residue, phi, psi]

        if fields:
            print fields
        print geometry

if __name__ == '__main__':
    # run as a program
    sys.exit(main(sys.argv))
