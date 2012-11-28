from gui import contactsDialog

def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
		'Visualize CASP RR contacts',
		label = 'CASP RR contacts viewer',
		command = lambda s=self : contactsDialog(s.root))
