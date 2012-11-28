from pymol import cmd
import sys, re, random, math, collections, operator

from util import *
from constants import *
from parser import *

def show_contacts(contactFile, target, chain, num_contacts, min_separation, contact_atom_mapping):
	"""Visualize contacts on a target"""


	with open(contactFile) as f_in:
		casprr = parse_casp_rr(f_in)

	sequence = get_sequence(target, chain) 
	contacts = casprr['contacts']

	print("casp: {0}\npdb:  {1}".format(sequence, casprr['sequence']))

	contacts = [ c for c in contacts if abs( c['i'] - c['j'] ) > min_separation]
	contacts = sorted(contacts, key=operator.itemgetter("conf"), reverse=True)
	contacts = contacts[0:num_contacts]


	constraints = generate_constraints(contacts, sequence, contact_atom_mapping)

	# clean up previously created objects
	if "contacts" in cmd.get_names(): cmd.delete("contacts")

	geom = []

	for constraint in constraints:

		# decorate dict with additional information for easy string formatting
		constraint['target'] = target
		constraint['chain'] = chain
	
		posx = atom_pos("/{target}//{chain}/{i}/{x}".format(**constraint))[0]
		posy = atom_pos("/{target}//{chain}/{j}/{y}".format(**constraint))[0]

		dst = math.sqrt((posx[0] - posy[0])**2 + (posx[1] - posy[1])**2 + (posx[2] - posy[2])**2)

		color = gradient_interpolate(dst)

		geom.extend( cylinder(posx, posy, c1=color, c2=color, r=CONSTRAINT_RADIUS) )

	cmd.load_cgo(geom, "contacts", 1)


#if __name__ == "__main__":
#	root = Tk()
#	contactsDialog(root)
#	root.mainloop()
