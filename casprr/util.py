from pymol import cmd, cgo
import math, collections, operator

from constants import *


def get_sequence(target, chain=""):
	"""Determine protein sequence from PDB ATOM fields"""

	if chain == "": chain = " "

	pdb = cmd.get_pdbstr(target).split("\n")

	seq3 = [line[17:20] for line in pdb if line[13:16] == "CA " and line[21] == chain]
	seq1 = [three_to_one[code] for code in seq3 if code] 

	return "".join(seq1)



def gradient_interpolate(pos, breaks={8: (0,1,0), 12:(1,1,0), 23:(1,0,0)}) :

	pos = float(pos)

	brks = list(breaks.items())
	brks = sorted(brks, key=lambda b: b[0])
	
	brks.insert(0, (float('-inf'), brks[0][1]))
	brks.append((float('inf'), brks[-1][1]))

	i = next( i for i, b in enumerate(brks) if pos < b[0])

	bfr = brks[i-1]
	bto = brks[i]

	if all( f == t for f,t in zip(bfr[1], bto[1])): return bfr[1]

	k = float(pos - bfr[0]) / (bto[0] - bfr[0])

	return tuple( (1-k) * f + k*t for f,t in zip(bfr[1], bto[1]) )
	

def cylinder(x1, x2, c1=(0,1,0), c2=(0,1,0), r=0.5):
	"""Create a CGO cylinder
	
	:param x1: a (x,y,z) tuple giving the first coordinate in angstroms
	:param x2: a (x,y,z) tuple giving the second coordinate in angstroms
	:param c1: a (r,g,b) tuple giving the color of the first coordinate (range [0:1])
	:param c2: a (r,g,b) tuple giving the color of the second coordinate (range [0:1])
	:param r: the radius of the cylinder in angstroms

	"""
	if isinstance(x1, str):
		x1 = atom_pos(x1)
	
	if isinstance(x2, str):
		x2 = atom_pos(x2)

	return [cgo.CYLINDER, x1[0], x1[1], x1[2], x2[0], x2[1], x2[2], r, c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]]

def atom_pos(selection):
	return cmd.get_model(selection, 1).get_coord_list()

