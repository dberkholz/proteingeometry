#!/usr/bin/env python

import sys
from pgd.angles import *

file='test.res'

def main():
	dbdict = create_all_databases(databases)
	for line in open(file):
		words = line.split()
		residue = words[0]
		phi = int(float(words[1]))
		psi = int(float(words[2]))
		fields, geometry = get_geometry(dbdict, residue, phi, psi)
		if fields:
			print fields
		print geometry

if __name__ == '__main__':
	sys.exit(main())
