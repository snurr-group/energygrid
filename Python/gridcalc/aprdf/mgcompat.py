# Compatibility layer between my cifparse.py and Arun's use of Pymatgen
from cifparse import CifData, CifBox, CifAtoms
from cifparse import clean_uc, clean_float
import numpy as np

def repeat(array, times=1):
	# Repeats elements of a list a specified number of times in the same order as the list
	# repeat([1,2,3], 3) yields [1,1,1,2,2,2,3,3,3]
	out_array = []
	for element in array:
		out_array.extend(times * [element])
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

		dx = transforms[0].repeat(num_atoms)
		dy = transforms[1].repeat(num_atoms)
		dz = transforms[2].repeat(num_atoms)

		# Apply changes to copies of atoms, and scale UC appropriately
		new_data['_atom_site_fract_x'] = (new_data['_atom_site_fract_x'] + dx) / mx
		new_data['_atom_site_fract_y'] = (new_data['_atom_site_fract_x'] + dy) / my
		new_data['_atom_site_fract_z'] = (new_data['_atom_site_fract_z'] + dz) / mz
		new_data['_cell_length_a'] = clean_float(self.data['_cell_length_a']) * mx
		new_data['_cell_length_b'] = clean_float(self.data['_cell_length_b']) * my
		new_data['_cell_length_c'] = clean_float(self.data['_cell_length_c']) * mz
		for key in ['_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma', '_symmetry_space_group_name_H-M']:
			new_data[key] = self.data[key]

		# Apply changes back to the class, which will require reinitialization
		self.data = new_data
		self.box = CifBox(self.data)
		self.atoms = CifAtoms(self.data)
		self.set_aliases()

