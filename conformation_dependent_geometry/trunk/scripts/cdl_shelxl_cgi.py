#!/usr/bin/python
"""This script accepts a PDB file from the user, runs cdl-shelxl
to produce a list of CDL restraints for the PDB file, and returns
those restraints to the user in a fresh browser page.
"""
import cgi
import cgitb; cgitb.enable()
import os, sys

# Write the top of the output page.

print "Content-type: text/html"
print 
print "<HTML>"

form_field = "PDBfile"
form = cgi.FieldStorage()

if form_field not in form:
    print "<HEAD><TITLE> cdl-shelxl </TITLE></HEAD>"
    print "<BODY><H1>Error</H1>"
    print "Please fill select a PDB file. </BODY></HTML>"
    sys.exit(1)

if not form.has_key(form_field): 
    print "<HEAD><TITLE> cdl-shelxl </TITLE></HEAD>"
    print "<BODY><H1>Error</H1>"
    print "I don\'t know how this happened but the sending form does not contain the \'", form_field, "\' key. </BODY></HTML>"
    sys.exit(1)

fileitem = form[form_field]
if not fileitem.file: 
    print "<HEAD><TITLE> cdl-shelxl </TITLE></HEAD>"
    print "<BODY><H1>Error</H1>"
    print "No PDB file was uploaded to the server. </BODY></HTML>"
    sys.exit(1)

# form.getfirst(form_field) is the contents of the uploaded file.

print "<HEAD><TITLE> ", fileitem.filename, "</TITLE></HEAD>"
print "<BODY><H1> ", fileitem.filename, "</H1>" 

# Send the up-loaded file to the script for processing.

import cdl_shelxl

try:
    cdl_shelxl.process(fileitem.file, 'html')
finally:
    print "</BODY></HTML>"

sys.exit(0)
