from distutils.core import setup
setup(name='pgd',
      version='0.1',
      description='Package for accessing Protein Geometry Database info',
      author='Donnie Berkholz',
      author_email='berkhold@science.oregonstate.edu',
      packages=['pgd'],
      scripts=['scripts/pgd-angles'],
      package_data={'pgd': ['data/*.txt'] }
      )
