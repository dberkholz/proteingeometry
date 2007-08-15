from distutils.core import setup
name='pgd-utils'
version='0.2'
setup(name=name,
      version=version,
      description='Package for accessing Protein Geometry Database info',
      author='Donnie Berkholz',
      author_email='berkhold@science.oregonstate.edu',
      packages=['pgd'],
      scripts=['scripts/pgd-angles'],
      package_data={'pgd': ['data/*.txt'] },
      data_files=[('share/doc/' + name + '-' + version,
                   ['karplus-definitions.jpg',
                    'examples/pgd-lookup.py',
                    'examples/test.res'])]
      )
