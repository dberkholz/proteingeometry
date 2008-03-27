import angles
#import get_geom_rmsd

angles_parser = angles.optparse_setup()
angles_options = []
angles.optlist, angles.args = angles_parser.parse_args(angles_options)

#get_geom_rmsd_parser = get_geom_rmsd.optparse_setup()
#get_geom_rmsd_options = []
#get_geom_rmsd.optlist, get_geom_rmsd.args = get_geom_rmsd_parser.parse_args(get_geom_rmsd_options)
