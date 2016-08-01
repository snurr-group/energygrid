# Compatibility layer between my cifparse.py and Arun's use of Pymatgen
from cifparse import CifData, CifBox, CifAtoms
from cifparse import clean_uc, clean_float
import numpy as np

# Needed for gcd:
import fractions
import copy


def gcd(list_of_ints):
	"""Calculate the greatest common denominator/factor of a list"""
	if len(list_of_ints) == 0:
		return None
	elif len(list_of_ints) == 1:
		return list_of_ints[0]
	values = copy.deepcopy(list_of_ints)
	# Start with one factor, which will have a gcf of itself
	gcf = values.pop()

	for coeff in values:
		# Calculate the gcf of the current gcf and next value
		gcf = fractions.gcd(coeff, gcf)
	return gcf

def _calc_empirical_formula(molec_formula):
	# Derived from aprdf:elements.py
	gcf = gcd(molec_formula.values())
	empirical_formula = dict.fromkeys(molec_formula)
	for element in molec_formula:
		# gcf will be a factor of each element (by definition), so integer division is fine
		empirical_formula[element] = molec_formula[element] / gcf
	return empirical_formula

def _formula_to_str(generic_formula_dict, sep = None, hide_ones = False):
	if sep is None:
		sep = ""
	elements = sorted(generic_formula_dict.keys())  # Pymatgen uses a different ordering, but this is otherwise consistent
	out = ""
	first_element = True
	for element in elements:
		if not first_element:  # use sep to delimit element entries
			out += sep
		out += element
		count = generic_formula_dict[element]
		if not(hide_ones and count == 1):
			out += str(count)
		first_element = False
	return out

def repeat(array, times=1):
	# Repeats elements of a list a specified number of times in the same order as the list
	# repeat([1,2,3], 3) yields [1,1,1,2,2,2,3,3,3]
	out_array = []
	for element in array:
		#out_array.extend(times * [element])
		out_array.extend([element for x in xrange(times)])
	return out_array


class CifParser(CifData):
	# Derived class from my cifparse.py utility
	# Aims to provide a basic compatibility layer for the pymatgen CIF reading used by Arun's energy grid code
	# See Structure class at http://pymatgen.org/_modules/pymatgen/core/structure.html
	def __init__(self, filename):
		super(CifParser, self).__init__(filename)
		self.set_aliases()

	def set_aliases(self):
		self.lattice = self.box
		self.frac_coords = self.atoms.fpos  # pymatgen also returns an Nx3 numpy array
		self.num_sites = self.atoms.count
		self.species = self.atoms.element
		self._calculate_formulas()

	def get_structures(self):
		# WARNING: The base cifparse.py class by default only reads the first structure in the CIF file
		# Therefore, this function is just a stub so that .get_structures()[0] returns the expected value
		return [self]

	def make_supercell(self, cell_array):
		# Make a MOF supercell based on the integer multiplications provided

		# Get the counts before/after for bookkeeping
		mx = int(cell_array[0])
		my = int(cell_array[1])
		mz = int(cell_array[2])
		assert mx >= 1
		assert my >= 1
		assert mz >= 1
		num_atoms = self.atoms.count
		num_mult = mx * my * mz
		num_new = self.atoms.count * num_mult

		# Preallocate the new data
		new_data = dict()
		labels = ['_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z']
		for key in labels:
			# Make an exact copy repeated, e.g. [1,2,3] 3 times is [1,1,1,2,2,2,3,3,3]
			new_data[key] = clean_uc(self.data[key]).repeat(num_mult)
		for key in ['_atom_site_type_symbol', '_atom_site_label']:  # Can't convert strings to floats...
			new_data[key] = repeat(self.data[key], num_mult)
		if '_atom_site_occupancy' in self.data.keys():
			del self.data['_atom_site_occupancy']

		# Calculate the unit cells we need, based on cifparse.permutations
		transforms = np.zeros((num_mult, 3), int)
		index = 0
		for x in range(mx):  # have to add one to include the upper bound
			for y in range(my):
				for z in range(mz):
					transforms[index] = [x, y, z]
					index += 1
		transforms = transforms.T

		# We want to do the transforms n times for each atom, e.g.
		# [0, 1, 2] -> [0, 1, 2, 0, 1, 2]
		# This is the opposite convention as the `repeat` command
		dx = np.tile(transforms[0], num_atoms)
		dy = np.tile(transforms[1], num_atoms)
		dz = np.tile(transforms[2], num_atoms)

		# Apply changes to copies of atoms, and scale UC appropriately
		new_data['_atom_site_fract_x'] = (new_data['_atom_site_fract_x'] + dx) / mx
		new_data['_atom_site_fract_y'] = (new_data['_atom_site_fract_y'] + dy) / my
		new_data['_atom_site_fract_z'] = (new_data['_atom_site_fract_z'] + dz) / mz
		new_data['_cell_length_a'] = clean_float(self.data['_cell_length_a']) * mx
		new_data['_cell_length_b'] = clean_float(self.data['_cell_length_b']) * my
		new_data['_cell_length_c'] = clean_float(self.data['_cell_length_c']) * mz
		for key in ['_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma', '_symmetry_space_group_name_H-M']:
			new_data[key] = self.data[key]

		# Apply changes back to the class, which will require reinitialization
		if False:  # Old debugging commands to figure out the transform
			import sys
			np.set_printoptions(threshold=sys.maxint)
			print cell_array
			print transforms[0], transforms[1], transforms[2]
			print new_data
		self.data = new_data
		self.box = CifBox(self.data)
		self.atoms = CifAtoms(self.data)
		self.set_aliases()

	def _calculate_formulas(self):
		# Calculates the molecular and empirical formulas for the CIF
		# Derived from aprdf:elements.py
		self.molecular_formula = dict()
		for item in self.atoms.element:
			if item in self.molecular_formula:
				self.molecular_formula[item] += 1
			else:
				self.molecular_formula[item] = 1
		# Check that all atoms are accounted for.  Print out mismatch if this assertion fails
		assert self.atoms.count == sum(self.molecular_formula.values()), str(self.atoms.count) + " != " + str(
			sum(self.molecular_formula.values()))

		# Get list of atom types
		self.atom_types = self.molecular_formula.keys()
		# Calculate the empirical formula (simplest integer representation of molecular formula)
		self.empirical_formula = _calc_empirical_formula(self.molecular_formula)
		return self

	def __str__(self):
		# Returns an imitation of Pymatgen's string representation.
		# Used for Structure_Information.dat, etc.
		summary = "Structure Summary "
		summary += "(" + _formula_to_str(self.molecular_formula, " ", False) + ")\n"
		summary += "Reduced Formula: " + _formula_to_str(self.empirical_formula, "", True) + "\n"

		summary += "abc   :  %9.6f %9.6f %9.6f\n" % (self.lattice.a, self.lattice.b, self.lattice.c)
		summary += "angles:  %9.6f %9.6f %9.6f\n" % (self.lattice.alpha, self.lattice.beta, self.lattice.gamma)

		summary += "Sites (" + str(self.atoms.count) + ")\n"
		for site in xrange(self.atoms.count):
			summary += str(site + 1) + " "
			summary += self.atoms.element[site] + 5 * " "
			summary += "%.6f" % self.atoms.fx[site] + 5 * " "
			summary += "%.6f" % self.atoms.fy[site] + 5 * " "
			summary += "%.6f" % self.atoms.fz[site] + "\n"

		return summary



