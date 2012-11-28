from Tkinter import *
import tkFileDialog, tkMessageBox
from pymol import cmd

from parser import contact_atoms 
from casprr import show_contacts

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
	lstObject = Listbox(frmMain, selectmode=SINGLE, exportselection=0, height=6)


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
	catms = sorted(contact_atoms.keys())
	lstAtomMappings = Listbox(frmMain, selectmode=MULTIPLE, exportselection=0, height=4)
	for n in catms:
		lstAtomMappings.insert(END, n)
	lstAtomMappings.selection_set(1)
	lstAtomMappings.grid(row=4, column=1, columnspan=2, sticky=N+E+S+W)

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
			return False

		if not lstAtomMappings.curselection():
			tkMessageBox.showwarning("No Atom Selection", "Please specify at least one set of atoms to map contacts on!")
			return False
		
		return True

	def confirm():

		if not validate(): return 

		contactFile = vTarget.get()
		target, chain = objects[int(lstObject.curselection()[0])]
		num_contacts = vContacts.get()
		min_separation = vSeparation.get()

		#atom_mapping = contact_atoms[ vAtomMappings.get() ]

		atom_mapping = [ contact_atoms[catms[int(k)]] for k in lstAtomMappings.curselection() ]

		print atom_mapping

		show_contacts(contactFile, target, chain, num_contacts, min_separation, atom_mapping)


	btnOK = Button(frmButtons, text="Show", command=confirm)
	btnOK.pack(side=RIGHT)

	browseTarget()
