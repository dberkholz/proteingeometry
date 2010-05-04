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

import angles
#import get_geom_rmsd

angles_parser = angles.optparse_setup()
angles_options = []
angles.optlist, angles.args = angles_parser.parse_args(angles_options)

from convert_kernel_regressions import convert_dunbrack_database

# Make sure our PGD-formatted databases are always up-to-date
convert_dunbrack_database()

#get_geom_rmsd_parser = get_geom_rmsd.optparse_setup()
#get_geom_rmsd_options = []
#get_geom_rmsd.optlist, get_geom_rmsd.args = get_geom_rmsd_parser.parse_args(get_geom_rmsd_options)
