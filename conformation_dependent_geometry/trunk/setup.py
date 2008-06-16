name='pgd-utils'
version='0.4.1'

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
      packages=['pgd'],
      scripts=['scripts/pgd-angles',
               'pgd/get_pdb_rmsd_stats.sh'],
      package_data={'pgd': ['data/*.txt'] },
      data_files=[('share/doc/' + name + '-' + version,
                   ['karplus-definitions.jpg',
                    'examples/pgd-lookup.py',
                    'examples/test.res'])],
      cmdclass={'sdist': sdist},
      )
