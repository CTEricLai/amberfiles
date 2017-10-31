#!/usr/bin/python
import MD_IO

# First, read topology file
# ex: AHEYKA.top is a topology file 
A = MD_IO.resinfo(top="AHEYKA.top")

# MD_IO is written in OOP, so use OOP approach to retrive the information
# For example,
# print out summary of residue 1
A.summary(1)

# print out bond terms for residue 1
for a in A.bondterms(1):
    print a

# print out angle terms for residue 2
for a in A.angleterms(2):
    print a

# print out dihedral terms for residue 3
for a in A.dihedralterms(3):
    print a
