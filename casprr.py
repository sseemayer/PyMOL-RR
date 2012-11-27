from Tkinter import *
import tkFileDialog
import tkMessageBox
from pymol import cmd
from operator import itemgetter
import collections
import sys, urllib2, zlib, re
import random

contact_atoms_C = { "A": "CB", "C": "CB", "D": "CB", "E": "CB", "F": "CB", "G": "CA", "H": "CB", "I": "CB", "K": "CB", "L": "CB", "M": "CB", "N": "CB", "P": "CB", "Q": "CB", "R": "CB", "S": "CB", "T": "CB", "V": "CB", "W": "CB", "Y": "CB" }
contact_atoms_CA = { "A": "CA", "C": "CA", "D": "CA", "E": "CA", "F": "CA", "G": "CA", "H": "CA", "I": "CA", "K": "CA", "L": "CA", "M": "CA", "N": "CA", "P": "CA", "Q": "CA", "R": "CA", "S": "CA", "T": "CA", "V": "CA", "W": "CA", "Y": "CA" }
contact_atoms_functional = { "A": "CB", "C": "SG", "D": "OD1", "E": "OE1", "F": "CZ", "G": "CA", "H": "CE1", "I": "CD1", "K": "NZ", "L": "CD1", "M": "CE", "N": "OD1", "P": "CG", "Q": "OE1", "R": "NH1", "S": "OG", "T": "OG1", "V": "CG1", "W": "CH2", "Y": "OH" }

three_to_one = { "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O"}

distance_groups_default = collections.OrderedDict([(0, {"name": "contact_close", "color": "green"}), (8, {"name": "contact_proximal", "color": "yellow"}), (13, {"name": "contact_distant", "color": "red"})])


def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
		'Visualize CASP RR contacts',
		label = 'CASP RR contacts viewer',
		command = lambda s=self : contactsDialog(s.root))


def get_sequence(target, chain=""):
	"""Determine protein sequence from PDB ATOM fields"""

	if chain == "": chain = " "

	pdb = cmd.get_pdbstr(target).split("\n")

	seq3 = [line[17:20] for line in pdb if line[13:16] == "CA " and line[21] == chain]
	seq1 = [three_to_one[code] for code in seq3 if code] 

	return "".join(seq1)

def parse_casp_rr(contactFile):
	"""Parse CASP RR file into a list of dicts
	
	Every list item will represent one contact prediction in the form of a dict with the following indices:

		i	The index of the first residue
		j	The index of the second residue
		dmin	The minimal predicted distance between i and j
		dmax	The maximal predicted distance between i and j
		conf	The confidence of the above prediction
	"""

	contacts = []

	with open(contactFile) as f_in:
		imax = 0
		for line in f_in:
			l = re.split(r"\s+", line.strip())
			if l[0] in ["PFRMAT", "TARGET", "AUTHOR", "REMARK", "METHOD", "MODEL", "END"] or len(l) < 2:
				continue

			i, j, dmin, dmax, conf = int(l[0]), int(l[1]), float(l[2]), float(l[3]), float(l[4])

			contacts.append({
				 'i': i	
				,'j': j
				,'dmin': dmin
				,'dmax': dmax
				,'conf': conf
			})

	return contacts


def show_contacts(contactFile, target, chain="", num_contacts=50, min_separation=4, contact_atom_mapping=contact_atoms_functional, distance_groups=distance_groups_default):
	"""Visualize contacts on a target"""

	sequence = get_sequence(target, chain) 

	contacts = parse_casp_rr(contactFile)
	contacts = [ c for c in contacts if abs( c['i'] - c['j'] ) > min_separation]
	contacts = sorted(contacts, key=itemgetter("conf"), reverse=True)
	contacts = contacts[0:num_contacts]

	# clean up previously created objects
	for dgroup in ( v['name'] for v in distance_groups.values()):
		if dgroup in cmd.get_names(): cmd.delete(dgroup)

	for contact in contacts:

		# decorate dict with additional information for easy string formatting
		contact['target'] = target
		contact['chain'] = chain
		contact['x'] = contact_atom_mapping[ sequence[ contact['i'] - 1 ] ]
		contact['y'] = contact_atom_mapping[ sequence[ contact['j'] - 1 ] ]
	
		# compute distance - this creates a dtemp object as a side effect which we will delete later
		dst = cmd.distance("dtemp", "/{target}//{chain}/{i}/{x}".format(**contact), "/{target}//{chain}/{j}/{y}".format(**contact))

		# find appropriate distance group
		dgroup = [ v for k,v in distance_groups.items() if dst >= k ][-1]

		# draw line in appropriate distance group
		cmd.distance(dgroup['name'], "/{target}//{chain}/{i}/{x}".format(**contact), "/{target}//{chain}/{j}/{y}".format(**contact))


	# color distance groups, hide labels
	for k,v in distance_groups.items():
		if v['name'] in cmd.get_names(): 
			cmd.color(v['color'], v['name'])
			cmd.hide('labels', v['name'])

	# remove temp distance group
	cmd.delete("dtemp")


def contactsDialog(root):
	"""Create GUI"""


	PADDING=5

	win = Toplevel(root, width=400, height=600, padx=PADDING, pady=PADDING)
	win.resizable(0,0)
	win.title("Visualize CASP RR contacts")

	#### MAIN CONTROLS ####

	frmMain = Frame(win)
	frmMain.pack(fill=BOTH, expand=1)
	frmMain.columnconfigure(0, weight=2, pad=PADDING)
	frmMain.columnconfigure(1, weight=5, pad=PADDING)
	frmMain.columnconfigure(2, weight=1, pad=PADDING)
	frmMain.rowconfigure(2, weight=4)

	Label(frmMain, text="RR file:").grid(row=0,column=0,sticky=N+E)
	vTarget = StringVar(win)
	txtTarget = Entry(frmMain, textvariable=vTarget)
	txtTarget.grid(row=0, column=1, sticky=W+E)
	
	def browseTarget():
		cf = tkFileDialog.askopenfilename(parent=win, filetypes=[("CASP RR files", ".CASPRR")])
		if cf:
			vTarget.set(cf)

	cmdTarget = Button(frmMain, text="...", command=browseTarget)
	cmdTarget.grid(row=0, column=2, sticky=W+E)
	
	Label(frmMain, text="Target:").grid(row=1,column=0,sticky=N+E)
	lstObject = Listbox(frmMain, selectmode=SINGLE)


	objects = []
	molecules = ( n for n in cmd.get_names() if cmd.get_type(n) == "object:molecule")
	for n in molecules:
		for ch in cmd.get_chains(n):
			objects.append((n, ch))
			lstObject.insert(END, "{0}/{1}".format(n,ch))

	if objects: lstObject.selection_set(0)

	lstObject.grid(row=1, column=1, columnspan=2, sticky=N+E+S+W)

	Label(frmMain, text="Min. separation:").grid(row=2,column=0,sticky=E)
	vSeparation = IntVar(win)
	vSeparation.set(23)
	sclSeparation = Scale(frmMain, from_=0, to=100, orient=HORIZONTAL, variable=vSeparation)
	sclSeparation.grid(row=2, column=1, columnspan=2, sticky=W+E)

	Label(frmMain, text="Num. contacts:").grid(row=3,column=0,sticky=E)
	vContacts = IntVar(win)
	vContacts.set(25)
	sclContacts = Scale(frmMain, from_=1, to=500, orient=HORIZONTAL, variable=vContacts)
	sclContacts.grid(row=3, column=1, columnspan=2, sticky=W+E)


	Label(frmMain, text="Use atoms:").grid(row=4,column=0,sticky=E)
	atom_mappings = collections.OrderedDict([
		 ("C-alpha", contact_atoms_CA)
		,("C-beta/C-alpha", contact_atoms_C)
		,("Functional", contact_atoms_functional)
	])
	vAtomMappings = StringVar(win)
	vAtomMappings.set("C-beta/C-alpha")
	lstAtomMappings = OptionMenu(frmMain, vAtomMappings, *atom_mappings.keys())
	lstAtomMappings.grid(row=4, column=1, columnspan=2, sticky=W+E)

	#### BUTTONS ROW ####

	frmButtons = Frame(win)
	frmButtons.pack(fill=X, expand=1)

	btnCancel = Button(frmButtons, text="Close", command=lambda: win.destroy())
	btnCancel.pack(side=RIGHT, pady=PADDING)

	def validate():
		if not vTarget.get():
			tkMessageBox.showwarning("No CASP RR file", "Please specify a valid CASP RR file to visualize!")
			return False

		if not lstObject.curselection():
			tkMessageBox.showwarning("No Mapping Target", "Please specify a molecule to map contacts on!")
		
		return True

	def confirm():

		if not validate(): return 

		contactFile = vTarget.get()
		target, chain = objects[int(lstObject.curselection()[0])]
		num_contacts = vContacts.get()
		min_separation = vSeparation.get()

		atom_mapping = atom_mappings[ vAtomMappings.get() ]

		show_contacts(contactFile, target, chain, num_contacts, min_separation, atom_mapping)


	btnOK = Button(frmButtons, text="Show", command=confirm)
	btnOK.pack(side=RIGHT)

	browseTarget()

#if __name__ == "__main__":
#	root = Tk()
#	contactsDialog(root)
#	root.mainloop()
