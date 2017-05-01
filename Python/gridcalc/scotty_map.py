#Python script to calculate the fractional attractive zone on a bunch of mofs, report that 
# number and also output the grid energy values for every mof in a set

# This is a modified code where the grid is geneerated only on the unit cell
#rest all are energy wise the same hence the energy values can basically be 
# replicated like the framework itself

write_xyz = False
write_vtk = False

if write_vtk:
	from pyevtk.hl import gridToVTK  # To output vtk data -- Advanced module need to be installed

try:
	from pymatgen.io.cifio import CifParser  # try to use Pymatgen, if installed
except ImportError:
	from aprdf.mgcompat import CifParser  # Alternative to pymatgen for importing CIFs
# Needs PyCIFRW library installed.  Might also need to run `conda install mingw` in Windows

import numpy as np  # Numerical calculations -- Basic module
import os           # System operations -- Basic module
import glob         # Wildcard expansion -- Basic module
import imp
import datetime
import time
import warnings
import gzip  # Consider writing the metadata as one huge table, with the file gzip.open('data.txt.gz', 'wb')

if write_xyz:
	xyz_mod = imp.load_source('xyz','xyz.py')  # A single file to make xyz file writing easy

# Define LJ Potential
def lj(eps,sig,rsq):
	E = (4*eps) * (((sig)**12/(rsq)**6) - ((sig)**6/(rsq)**3)) 
	return E

# Helper function for making new directories
def mkdir_if_new(path):
	if not os.path.isdir(path):
		os.mkdir(path)

#-------------------------------------------------------------------------------------------------------------
# Read the forcefield information from the RASPA style forcefield definition
#-------------------------------------------------------------------------------------------------------------

# Define the probe LJ parameters (H2 for the time being--change here)
eps1 = 37.3  # in kelvin
sig1 = 3.31  # in angstrom
rcut = 12.8  # Maximum range under consideration for the LJ interactions

# Units conversion from Kelvin to J/mol (to use integer values)
K_TO_J_MOL = 8.314

# Set the cut offs for defining the energy metric
e_low=-7000 # in kelvin units
e_high=-116

# Grid spacing
grid_spacing = 1.0  # One grid point per angstrom approximately

# Histogram parameters in [J/mol].  bin_min is automatically calculated from the minimum in the dataset
bin_max = 5000
bin_width = 10

# Read the forcefield information from RASPA force field directory
# I made some changes there to make life easy
force_field = 'UFF'
#pa_file = open('./forcefield/'+forcefield+'/pseudo_atoms.def').readlines()
#NumberOfPseudoAtoms = int(pa_file.[1].split()[0])

ff_file = open('./forcefield/'+force_field+'/force_field_mixing_rules.def').readlines()

# --------------------------------------------------------------------------------------------------------------------
# Prepare dir structure
# ---------------------------------------------
path_orig = os.getcwd()
os.chdir('CIF_FILES')
path_files=os.getcwd()
cif_list=glob.glob('*.cif')  # Assumes lowercase .cif, like cif_file_name below
cif_list = [x[:-4] for x in cif_list]  # Strip off the suffix

#-------------------------------------------------------------------------------------------------------------
# The heart of the code
#-------------------------------------------------------------------------------------------------------------
# Set up output files
mkdir_if_new('Grids/')
mkdir_if_new('Histograms/')
mkdir_if_new('Stats/')
details_file=open('Stats/Details.txt','w')
details_file.write('cif\tnx_total\tny_total\tnz_total\tnx_cells\tny_cells\tnz_cells\tlx\tly\tlz\talpha\tbeta\tgamma\tvf?\n')
metric_summary=open('Stats/Metric.txt','w')
metric_summary.write('cif\tmetric\n')
timer_file=open('Stats/timer.txt','w')
timer_file.write('cif\ttime (s)\n')
missing_params=open('Stats/bad_ff.txt', 'w')
open_stat_files = [details_file, metric_summary, timer_file, missing_params]

# write timestamp
starttime = time.clock()
fout = open('Stats/timestamp.txt','w')
fout.write('Timestamp: {:%Y-%b-%d %H:%M:%S}\n'.format(datetime.datetime.now()))

for name_index in range(len(cif_list)):
	timer_start = time.time()
	cif_file_name = cif_list[name_index] + '.cif'
	print cif_file_name
	
	# Read the coordinates from the cif file using pymatgen
	f1=CifParser(cif_file_name)
	struct=f1.get_structures()[0]
	
	
	# Need the box lengths to ensure the distance criterion
	# Create the cell matrix to ensure this criterion
	
	lx_unit = struct.lattice.a
	ly_unit = struct.lattice.b
	lz_unit = struct.lattice.c
	
	alpha = struct.lattice.alpha * (np.pi/180.0)
	beta  = struct.lattice.beta * (np.pi/180.0)
	gamma = struct.lattice.gamma * (np.pi/180.0)
	
	
	
	unit_cell_matrix = [[1 , 0, 0],[np.cos(gamma),np.sin(gamma),0],[0, 0, 1]]
	unit_rec_space = np.linalg.inv(unit_cell_matrix)
	
	
	# Make super cell based on the cut off criteria
	# the minimum box edge must be greater than twice the maximum range of interaction  
	# set to fixed 14 angstroms
	
	nx_cells = np.ceil(2*rcut/lx_unit)    # magic formula
	ny_cells = np.ceil(2*rcut/ly_unit)
	nz_cells = np.ceil(2*rcut/lz_unit) 
	
	
	struct.make_supercell([nx_cells,ny_cells,nz_cells]) # Structure is made into a super cell
	coord = np.array(struct.frac_coords)  # The whole thing scaled to [0,1] in all D's
	Number_Of_Atoms = struct.num_sites
	
	
	# Redefine the box matrix since we made a supercell
	lx = struct.lattice.a
	ly = struct.lattice.b
	lz = struct.lattice.c
	
	alpha = struct.lattice.alpha * (np.pi/180.0)
	beta  = struct.lattice.beta * (np.pi/180.0)
	gamma = struct.lattice.gamma * (np.pi/180.0)
	
	A = [[lx , 0, 0],[ly*np.cos(gamma),ly*np.sin(gamma),0],[0, 0, lz]]
	A = np.array(A).T
	A_inv = np.linalg.inv(A)
	
	# Define grid points in a unit box (Non adaptive grid)
	# Number of grid points
	nx = int(lx_unit/grid_spacing)
	ny = int(ly_unit/grid_spacing)
	nz = int(lz_unit/grid_spacing)
	
	# Intially the grids are defined on a only on the unit cell, which is only a tiny part of the unit box
	# The unit box corresponds to the entire super cell
	x_grid = np.linspace(0,1.0/nx_cells,nx)
	y_grid = np.linspace(0,1.0/ny_cells,ny)
	z_grid = np.linspace(0,1.0/nz_cells,nz)
	
	
	
	# Read the corresponding forcefield parameters (pseudo atoms information)
	# for the atoms in the MOF
	
	#Get atom names in the MOF
	mof_atm_names = []
	for i in range(Number_Of_Atoms):
		mof_atm_names.append(str(struct.species[i]))
	
	
	# Load the LJ parameters from the force fields
	ff_sig = dict()
	ff_eps = dict()
	for line in ff_file:
		fields = line.strip().split()
		if len(fields) == 4 and fields[1] == "lennard-jones":
			ff_sig[fields[0]] = fields[3]
			ff_eps[fields[0]] = fields[2]

	# Ensure that all atom types are defined
	mof_atm_indices = np.zeros((Number_Of_Atoms, 1))
	for i, element in enumerate(mof_atm_names):
		if (not ff_sig.has_key(element)) or (not ff_eps.has_key(element)):
			warnings.warn("Atom type " + element + " not defined in forcefield", UserWarning)
			missing_params.write(cif_file_name + '\t' + element + '\n')
			ff_sig[element] = 0.0
			ff_eps[element] = 0.0


	pot=np.zeros((nx,ny,nz)) # this is over just the unit cell
	
	#Calculate LJ interaction energy
	for atm_index in range(Number_Of_Atoms):
		# Lorentz Berthelot
		sig2 = float(ff_sig[mof_atm_names[atm_index]])
		eps2 = float(ff_eps[mof_atm_names[atm_index]])

		sig = (sig1 + sig2) / 2
		eps = np.sqrt(eps1 * eps2)
		for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					grid_point = np.array([x_grid[i],y_grid[j],z_grid[k]]) # fractional grid coordinater unit box
					# Don't convert back to Cartesian coordinates until we apply PBC.  Work in fractional coords
					# grid_point = np.dot(A , grid_point) # Cartesian coordinates
					
					drij = grid_point - coord[atm_index]
					# Apply PBC
					drij -= drij.round()
					# Minor detail: python's round appears to round towards +/- inf; numpy is round to nearest even
					
					#Transform back to find the actual distance
					r_act = np.dot(A , drij)
					rsq=r_act[0]**2+r_act[1]**2+r_act[2]**2
					
					if np.sqrt(rsq) <= rcut:
						pot[i][j][k] += lj(eps,sig,rsq)
	
	
	# Unit conversion
	pot *= K_TO_J_MOL
	
	#--------------Replicate the potential array so as to mimic the super cell
	pot_repeat = np.tile(pot,(nx_cells,ny_cells,nz_cells))
	
	
	nx_total = int(nx*nx_cells)
	ny_total = int(ny*ny_cells)
	nz_total = int(nz*nz_cells)
	
	
	#-------------------------------------------------------------------------------------------------------------
	# Output the VTS, grid energy values and attractive zone
	#-------------------------------------------------------------------------------------------------------------
	
	N_grid_total = nx_total *ny_total *nz_total
	
	#Write the VTK file
	if write_vtk:
		# Define the crazy unstructured grid
		dx, dy, dz = 1.0/nx_total, 1.0/ny_total, 1.0/nz_total
		X = np.arange(0, 1 + 0.1*dx, dx, dtype='float64')
		Y = np.arange(0, 1 + 0.1*dy, dy, dtype='float64')
		Z = np.arange(0, 1 + 0.1*dz, dz, dtype='float64')

		x = np.zeros((nx_total , ny_total , nz_total))
		y = np.zeros((nx_total , ny_total , nz_total))
		z = np.zeros((nx_total , ny_total , nz_total))

		for k in range(nz_total):
			for j in range(ny_total):
				for i in range(nx_total):
					x[i,j,k] = X[i]
					y[i,j,k] = Y[j]
					z[i,j,k] = Z[k]

		for k in range(nz_total):
			for j in range(ny_total):
				for i in range(nx_total):
					[x[i,j,k], y[i,j,k],z[i,j,k]] = np.dot(A,[x[i,j,k], y[i,j,k],z[i,j,k]])

		# Write pot into .vts file
		gridToVTK("./pot", x, y, z, pointData = {"Potential" : pot_repeat})
	
	
	# Histogram into predefined bins
	N_grid_inner = nx * ny * nz
	e_vals = np.reshape(pot,(N_grid_inner,1))  # Reshape into linear array, without supercell duplication

	if min(e_vals) > 0:
		print 'Nonporous material: no attractive region. ', cif_file_name
		bin_min = 0
	else:
		bin_min = int(min(e_vals)/bin_width) * bin_width - bin_width  # Align the lowest bin with bin_width
	bins1 = range(bin_min, bin_max+bin_width, bin_width)
	bins1.extend([np.inf])
	e_hist, binedges1 = np.histogram(e_vals,bins=bins1) # Histogram
	bin_left = np.asarray(binedges1[:-1], 'int')
	bin_right = np.asarray(binedges1[1:], 'int')
	bin_right[-1] = bin_max * 10  # Pseudo final bin indicates all the energies above the set maximum
	data = np.transpose(np.vstack((bin_left, bin_right, e_hist)))
	np.savetxt('Histograms/' + cif_list[name_index] + '_histogram.txt.gz', data, fmt='%d')  # Write histogram data
	vf = float(sum(np.less(e_vals, 0))) / N_grid_inner  # Approximating void frac as points with E < 0


	# Write details about the unit cells and replications
	stuff = [str(x) for x in [cif_file_name, nx_total, ny_total, nz_total, nx_cells, ny_cells, nz_cells, lx, ly, lz, alpha, beta, gamma, vf]]
	details_file.write('\t'.join(stuff) + '\n')

	# Write the raw energy values.  Could also consider a 3d npy file, but 3D arrays are less generally compatible
	np.savetxt('Grids/'+cif_list[name_index]+'_Energy_Values.txt.gz', e_vals, fmt='%.6g')
	
	
	#Print the fraction of the attractive zone
	metric_summary.write(cif_file_name+'\t')
	metric_summary.write(str(float(sum((e_vals < e_high) & (e_vals> e_low))[0])/N_grid_inner) + '\n')
	
	# write the raw grid data into text file straight array
	# Write the corresponding xyz coordinates
	# Write the corresponding MOF coordinates to .xyz file using xyz.py
	if write_xyz:
		f4=open(cif_list[name_index]+'.xyz','w')
		coord = np.array(struct.frac_coords)
		out_coord = np.zeros((np.shape(coord)))
		for i in range(len(coord)):
			out_coord[i]=np.dot(A,coord[i])
		xyz_mod.write_xyz(f4,out_coord,title=cif_list[name_index]+'.xyz',atomtypes=mof_atm_names)
		f4.close()

	# Write an approximate elapsed time (seconds) from analyzing this CIF file
	timer_end = time.time()
	elapsed_seconds = timer_end - timer_start
	timer_file.write(cif_file_name + '\t' + str(elapsed_seconds) + '\n')

	# Write files to disk to make sure we don't lose any progress if the script is closed or crashes
	[f.flush() for f in open_stat_files]
	[os.fsync(f.fileno()) for f in open_stat_files]

# Clean up output files
[f.close() for f in open_stat_files]

#print timing
endtime = time.clock()
elaptime = endtime-starttime
fout.write('Timestamp: {:%Y-%b-%d %H:%M:%S}\n'.format(datetime.datetime.now()))
fout.write('Elapsed time\t%f\n' % (elaptime))
fout.close()

os.chdir(path_orig)


