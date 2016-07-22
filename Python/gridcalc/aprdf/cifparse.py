import math
import numpy as np
import CifFile
import re


def permutations(combo_array):
	# Get permutations for a 3 length array with each element going
	# from -n to +n.
	
	# Extract the coordinates
	xc = combo_array[0]
	yc = combo_array[1]
	zc = combo_array[2]
	
	# Each element will have 2n + 1 possibilities, including zero
	num_combos = (2*xc + 1) * (2*yc + 1) * (2*zc + 1)
	combos = np.zeros((num_combos, 3), int)
	
	index = 0
	for x in range(-1*xc, xc+1):  # have to add one to include the upper bound
		for y in range(-1*yc, yc+1):
			for z in range(-1*zc, zc+1):
				combos[index] = [x, y, z]
				index = index + 1
	
	return combos


# Section with some new utility functions for increased robustness
# TODO: also raise an Exception unless the unit cell is P1
def wrap_uc(dfpos):
	# Apply PBC except in (0,1) instead of (-0.5, 0.5)
	# Python modulo of negatives/floating point is not unexpected
	# e.g. see http://python-history.blogspot.com/2010/08/why-pythons-integer-division-floors.html
	# Utility function copied from Ben's RCSR conversion code
	dfr = dfpos % 1  # only want the fractional, positive component
	# e.g. np.array([0.5, 1.1, -0.2], float) % 1 is what we want
	return dfr

def clean_float(string_number):
	# Converts a string to a float, removing any parentheses from uncertainty
	if type(string_number) in [str, unicode]:
		string_number = re.sub(r'\([^\)]*\)', '', string_number)
	return float(string_number)

def clean_uc(array_with_uncertainty):
	# clean_float on an array, and ensure that coordinates are in [0, 1)
	return wrap_uc(np.array([clean_float(i) for i in array_with_uncertainty], float))



class CifBox:
	# Simulation box parameters
	def __init__(self, cifdata, crystal_radius=0):
		if cifdata['_symmetry_space_group_name_H-M'] not in ["P1", "P 1"]:
			raise TypeError("CIF parsing is currently only supported for P1 unit cells")
		# By default, use standard periodic boundary conditions on a generic triclinic cell
		# If crystal_radius is defined (and larger than one Angstrom), use ghost atoms
		self.a = clean_float(cifdata['_cell_length_a'])  # read number instead of unicode string
		self.b = clean_float(cifdata['_cell_length_b'])
		self.c = clean_float(cifdata['_cell_length_c'])
		self.alpha = clean_float(cifdata['_cell_angle_alpha'])
		self.beta = clean_float(cifdata['_cell_angle_beta'])
		self.gamma = clean_float(cifdata['_cell_angle_gamma'])
		self.calculate_cell(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
		if crystal_radius > 1.0:
			# Precompute and automatically choose the proper pbc_cart_dists function based on circumstance
			self._precompute_symmetry_ops(crystal_radius)
			self.pbc_cart_dists = self._cart_crystal_dist
		else:
			self.pbc_cart_dists = self._standard_pbc_cart_dist
	
	def calculate_cell(self, a, b, c, alpha, beta, gamma):
		# Calculate unit cell matrix for triclinic coordinate conversions
		DEG2RAD = math.pi/180
		alpha = alpha * DEG2RAD
		beta = beta * DEG2RAD
		gamma = gamma * DEG2RAD
		cos = math.cos
		sin = math.sin
		
		# Volume of the unit cell
		self.volume = a*b*c * math.sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*cos(alpha)*cos(beta)*cos(gamma))
		
		# projections from fractional coords to Cartesian
		self.ax = a
		self.ay = 0.0
		self.az = 0.0
		self.bx = b * cos(gamma)
		self.by = b * sin(gamma)
		self.bz = 0.0
		self.cx = c * cos(beta)
		self.cy = c * (cos(alpha) - cos(gamma)*cos(beta)) / sin(gamma)
		self.cz = math.sqrt(c**2 - self.cx**2 - self.cy**2)  # distance formula
		
		# Transpose the transform array so that the matrix multiplication fract_coords*cell_transf=cart_coords
		# Could transpose with .T if we wanted multiplication in the other order
		self.cell_transform = np.array([[self.ax, self.ay, self.az], [self.bx, self.by, self.bz], [self.cx, self.cy, self.cz]], float)
		# Get the inverse transform for Cartesian to fractional coordinates
		self.inv_transform = np.linalg.inv(self.cell_transform)
		return self
	
	def apply_pbc(self, dfpos):
		# Apply periodic boundary conditions to fractional coordinate or delta
		# Input and output are position/displacement vectors in fractional coordinates
		# dfx = dfpos[0]
		# dfy = dfpos[1]
		# dfz = dfpos[2]
		#  Apply minimum image convention
		# dfx = dfx - round(dfx)
		# dfy = dfy - round(dfy)
		# dfz = dfz - round(dfz)
		# dfr = np.array([dfx, dfy, dfz])
		# Rewrite in numpy idiom for about 40% speed boost
		dfr = dfpos - dfpos.round()
		# Minimum image convention (in most cases?) is same for orthorhomic
		# and triclinic systems
		return dfr
	
	def frac_to_cart(self, fpos):
		# Convert fractional coordinates/delta to Cartesian coordinates
		# x = self.ax * fpos[0] + self.bx * fpos[1] + self.cx * fpos[2]
		# y = self.ay * fpos[0] + self.by * fpos[1] + self.cy * fpos[2]
		# z = self.az * fpos[0] + self.bz * fpos[1] + self.cz * fpos[2]
		# cpos = np.array([x, y, z])
		# Rewrite using numpy dot product, 3x1 times 3x3 transform array, equivalent to above
		cpos = np.dot(fpos, self.cell_transform)
		return cpos
	
	def pbc_cart_dists(self, fpos):
		# Empty function that gets redefined by class initiation to one of the next two
		pass
	
	def _standard_pbc_cart_dist(self, dfpos):
		dr = self.frac_to_cart(self.apply_pbc(dfpos))
		rij = math.sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])
		return np.array([rij], float)
	
	def _cart_crystal_dist(self, dfpos):
		# Calculates all potential lengths to ghost atoms within a length cutoff		
		# First convert dfpos to periodic boundary conditions.  Then look for ghosts.
		dfr = dfpos - dfpos.round()
		
		# Apply symmetry operations.  Do not apply the cell length since Cartesian conversion is next
		# Like standard PBC, symmetry should start by subtracting box lengths from positive
		# values and adding box lengths to negative values.
		# dfr_sign = np.array([1,1,1], int)
		# dfr_sign[dfr < 0] = -1
		# frac_lengths = dfr - dfr_sign * self.symmetry_ops
		# Since symmetry_ops is symmetric around zero, subtracting and adding are equivalent.  Skip computing the sign to save >10% computation time
		frac_lengths = dfr - self.symmetry_ops
		
		# Convert from fractional coordinates to Cartesian
		# For now, just apply the fractional to cartesian conversion with the general method
		drs = np.dot(frac_lengths, self.cell_transform)
		rijs = np.sqrt(np.sum(drs**2, axis=1))  # sum each row of squared terms
		lengths = rijs[rijs<=self.radius]
		return lengths
	
	def _precompute_symmetry_ops(self, radius):
		# Calculate the max number of unit cell images (brute force)
		# Based on worst case of a vector entirely in one direction of the unit cell.
		# symmetry_ops represent ghost atom(s) that may be necessary if the desired radius extends beyond
		# half a box length (which loops back due to minimum image convention)
		self.radius = radius
		cell_dims = np.array([self.a, self.b, self.c], float)
		# Add a small numerical buffer to make sure borderline 0.5 rounds up
		# Worst case, this only adds computation time, since extra-long distances are automatically filtered out later
		images = np.array((radius / cell_dims + 0.01).round(), int)
		self.symmetry_ops = permutations(images)
		self.num_symmetry_ops = len(self.symmetry_ops)


class CifAtoms:
	def __init__(self, cifdata):
		self.fx = clean_uc(cifdata['_atom_site_fract_x'])
		self.fy = clean_uc(cifdata['_atom_site_fract_y'])
		self.fz = clean_uc(cifdata['_atom_site_fract_z'])
		self.fpos = np.array([self.fx,self.fy,self.fz]).astype(np.float).T	# n x 3 matrix of positions
		self.element = cifdata['_atom_site_type_symbol']
		self.label = cifdata['_atom_site_label']
		self.count = len(self.element)
		# Read optional parameters
		self.occupancy = None
		if '_atom_site_occupancy' in cifdata.keys():
			self.occupancy = cifdata['_atom_site_occupancy']


class CifData(object):
	# Class that reads in a single block of data from a CIF file and organizes
	# the atom and simulation box variables
	# Makes initialization as easy as mydata = CifData(filename) and usage as
	# mydata.atoms.fpos[10][2] for the fractional z coordinate of the 11th atom
	def __init__(self, filename, crystal_radius=0):
		self.data = CifFile.ReadCif(filename).first_block()
		self.box = CifBox(self.data, crystal_radius)
		self.atoms = CifAtoms(self.data)
