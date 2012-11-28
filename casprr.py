from Tkinter import *
import tkFileDialog, tkMessageBox
from pymol import cmd, cgo
import sys, urllib2, zlib, re, random, math, collections, operator

contact_atoms_C = { "A": "CB", "C": "CB", "D": "CB", "E": "CB", "F": "CB", "G": "CA", "H": "CB", "I": "CB", "K": "CB", "L": "CB", "M": "CB", "N": "CB", "P": "CB", "Q": "CB", "R": "CB", "S": "CB", "T": "CB", "V": "CB", "W": "CB", "Y": "CB" }
contact_atoms_CA = { "A": "CA", "C": "CA", "D": "CA", "E": "CA", "F": "CA", "G": "CA", "H": "CA", "I": "CA", "K": "CA", "L": "CA", "M": "CA", "N": "CA", "P": "CA", "Q": "CA", "R": "CA", "S": "CA", "T": "CA", "V": "CA", "W": "CA", "Y": "CA" }
contact_atoms_functional = { "A": "CB", "C": "SG", "D": "OD1", "E": "OE1", "F": "CZ", "G": "CA", "H": "CE1", "I": "CD1", "K": "NZ", "L": "CD1", "M": "CE", "N": "OD1", "P": "CG", "Q": "OE1", "R": "NH1", "S": "OG", "T": "OG1", "V": "CG1", "W": "CH2", "Y": "OH" }

three_to_one = { "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O"}

CONSTRAINT_RADIUS = 0.1

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

def show_contacts(contactFile, target, chain="", num_contacts=50, min_separation=4, contact_atom_mapping=contact_atoms_functional):
	"""Visualize contacts on a target"""

	sequence = get_sequence(target, chain) 

	contacts = parse_casp_rr(contactFile)
	contacts = [ c for c in contacts if abs( c['i'] - c['j'] ) > min_separation]
	contacts = sorted(contacts, key=operator.itemgetter("conf"), reverse=True)
	contacts = contacts[0:num_contacts]

	# clean up previously created objects
	if "contact" in cmd.get_names(): cmd.delete("contact")

	geom = []

	for contact in contacts:

		# decorate dict with additional information for easy string formatting
		contact['target'] = target
		contact['chain'] = chain
		contact['x'] = contact_atom_mapping[ sequence[ contact['i'] - 1 ] ]
		contact['y'] = contact_atom_mapping[ sequence[ contact['j'] - 1 ] ]
	

		posx = atom_pos("/{target}//{chain}/{i}/{x}".format(**contact))[0]
		posy = atom_pos("/{target}//{chain}/{j}/{y}".format(**contact))[0]

		dst = math.sqrt((posx[0] - posy[0])**2 + (posx[1] - posy[1])**2 + (posx[2] - posy[2])**2)

		color = gradient_interpolate(dst)

		geom.extend( cylinder(posx, posy, c1=color, c2=color, r=CONSTRAINT_RADIUS) )

	cmd.load_cgo(geom, "contact", 1)


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
			print(n,ch)
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
