# -*- coding: utf-8 -*-
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

name='pgd-utils'
version='0.5'

import os
from distutils.core import setup
from distutils.command.sdist import sdist as _sdist

class sdist(_sdist):
    """Specialized source distribution creation"""
    if os.path.exists('.git'):
        cmd='git log > ChangeLog'
        os.system(cmd)

setup(name=name,
      version=version,
      description='Package for accessing Protein Geometry Database info',
      author='Donnie Berkholz',
      author_email='berkhold@science.oregonstate.edu',
      packages=['conformation_dependent_geometry'],
      scripts=['scripts/pgd-angles',
               'scripts/cdl-shelxl.py',
               'scripts/convert_pdb_to_refmac_restraints.py',
               'scripts/get_pdb_rmsd_stats.sh'],
      package_data={'pgd': ['data/*.bz2'] },
      data_files=[('share/doc/' + name + '-' + version,
                   ['karplus-definitions.jpg',
                    'examples/pgd-lookup.py',
                    'examples/test.res'])],
      cmdclass={'sdist': sdist},
      )
