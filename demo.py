#!/usr/bin/python

import MD_IO

A = MD_IO.resinfo(top='AHEYKA.top')

# summary of residue 1
A.summary(1)

# bond terms of residue 2
for a in A.bondterms(2):
    print a

# angle terms of residue 3
for a in A.angleterms(3):
    print a

# dihedral terms of residue 4
for a in A.dihedralterms(4):
    print a
