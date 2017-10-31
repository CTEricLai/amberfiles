#!/usr/bin/python

import MD_IO

A = MD_IO.resinfo(top="AHEYKA.top")

A.summary(2)

for a in A.dihedralterms(1):
	print a
