import re, collections, collections_extra

contact_atoms = {
	 "CA": { "A": "CA", "C": "CA", "D": "CA", "E": "CA", "F": "CA", "G": "CA", "H": "CA", "I": "CA", "K": "CA", "L": "CA", "M": "CA", "N": "CA", "P": "CA", "Q": "CA", "R": "CA", "S": "CA", "T": "CA", "V": "CA", "W": "CA", "Y": "CA" }
	,"CB/CA": { "A": "CB", "C": "CB", "D": "CB", "E": "CB", "F": "CB", "G": "CA", "H": "CB", "I": "CB", "K": "CB", "L": "CB", "M": "CB", "N": "CB", "P": "CB", "Q": "CB", "R": "CB", "S": "CB", "T": "CB", "V": "CB", "W": "CB", "Y": "CB" }
	,"functional": { "A": "CB", "C": "SG", "D": "OD1", "E": "OE1", "F": "CZ", "G": "CA", "H": "CE1", "I": "CD1", "K": "NZ", "L": "CD1", "M": "CE", "N": "OD1", "P": "CG", "Q": "OE1", "R": "NH1", "S": "OG", "T": "OG1", "V": "CG1", "W": "CH2", "Y": "OH" }
}

def parse_casp_rr(f_in):
	"""Parse CASP RR file into python dict.
	
	The dict will have the following elements:
	
	- target: The target name, e.g. T1234
	- author: The author id, eg. 1234-5678-9000
	- remark: A list of REMARK lines with the REMARK removed
	- method: A list of METHOD lines with the METHOD removed
	- sequence: The protein sequence of the target
	- contacts: A list of dicts. Every list item will represent one contact prediction in the form of a dict with the following indices:
		i	The index of the first residue
		j	The index of the second residue
		dmin	The minimal predicted distance between i and j
		dmax	The maximal predicted distance between i and j
		conf	The confidence of the above prediction

	The parser will do some state checking to ensure it parses a valid CASP RR file.
	
	"""

	RE_PFRMAT = re.compile(r'^PFRMAT\s+RR\s*$')
	RE_TARGET = re.compile(r'^TARGET\s+(.*?)\s*$')
	RE_AUTHOR = re.compile(r'^AUTHOR\s+(.*?)\s*$')
	RE_REMARK = re.compile(r'^REMARK\s+(.*?)\s*$')
	RE_METHOD = re.compile(r'^METHOD\s+(.*?)\s*$')
	RE_MODEL  = re.compile(r'^MODEL\s+(\d+)\s*$')
	RE_PRED   = re.compile(r'^\d')
	RE_SPLIT  = re.compile(r'\s+')
	RE_END  = re.compile(r'^END\s*$')

	COL_MODIFIERS = collections.OrderedDict([
		 ("i",    int)
		,("j",    int)
		,("dmin", float)
		,("dmax", float)
		,("conf", float)
	])

	target = None
	author = None
	remark = []
	method = []
	
	sequence = []
	contacts = []

	state = 0	# 0: header 1: sequence 2: contact 3: eof

	for i, line in enumerate(f_in):

		line = line.strip()
		if line == "": next

		if state == 0: # header
			if RE_PFRMAT.match(line): 
				pass
			elif RE_TARGET.match(line):
				target = RE_TARGET.match(line).group(1)
			elif RE_AUTHOR.match(line): 
				author = RE_AUTHOR.match(line).group(1)
			elif RE_REMARK.match(line):
				remark.append(RE_REMARK.match(line).group(1))
			elif RE_METHOD.match(line):
				method.append(RE_METHOD.match(line).group(1))
			elif RE_MODEL.match(line):
				state = 1
			else:
				raise Exception("Unrecognized header element in line {i}: '{line}'".format(i=i, line=line))

		elif state == 1: # sequence
			if not RE_PRED.match(line):
				sequence.append(line) 
			else:
				state = 2

		elif state == 2: # contact

			if not RE_END.match(line):

				l = RE_SPLIT.split(line)

				c = { key: COL_MODIFIERS[key](l[i]) for i, key in enumerate(COL_MODIFIERS.keys()) }	
				contacts.append(c)

			else:
				state = 3

		elif state == 3:
			raise Exception("Commands after END in line {i}: '{line}'".format(i=i, line=line))

	if state != 3:
		raise Exception("No END at EOF!")

	return {
		 'target': target
		,'author': author
		,'remark': remark
		,'method': method
		,'sequence': "".join(sequence)
		,'contacts': contacts
	}


def generate_constraints(contacts, sequence, constraint_atoms=[contact_atoms['CB/CA'], contact_atoms['functional'], contact_atoms['CA']], dist_override=None):
	
	constraints = []
	for c in contacts:

		for catm in constraint_atoms:

			constraint = c.copy()	

			distance = (constraint['dmin'] + constraint['dmax']) / 2
			dminus = dplus = distance - constraint['dmin']

			if dist_override: 
				distance, dminus, dplus = dist_override

			constraint.update({
				 "x": 		catm[sequence[c['i']-1]]
				,"y": 		catm[sequence[c['j']-1]]
				,"distance":	distance	
				,"dminus": 	dminus
				,"dplus":	dplus
			})
			
			constraints.append(constraint)

	# remove duplicate constraints
	constraints = [ dict(y) for y in collections_extra.OrderedSet( tuple(x.items()) for x in constraints )]

	return constraints
