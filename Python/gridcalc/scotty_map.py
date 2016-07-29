#Python script to calculate the fractional attractive zone on a bunch of mofs, report that 
# number and also output the grid energy values for every mof in a set

from pyevtk.hl import gridToVTK  # To output vtk data -- Advanced module need to be installed

try:
	from pymatgen.io.cifio import CifParser  # try to use Pymatgen, if installed
except ImportError:
	from aprdf.mgcompat import CifParser  # Alternative to pymatgen for importing CIFs
# Needs PyCIFRW library installed.  Might also need to run `conda install mingw` in Windows
import numpy as np  # Numerical calculations -- Basic module
import os           # System operations -- Basic module


# Define LJ Potential
def lj(eps,sig,rsq):
	E = (4*eps) * (((sig)**12/(rsq)**6) - ((sig)**6/(rsq)**3)) 
	return E



#-------------------------------------------------------------------------------------------------------------
# Read the forcefield information from the RASPA style forcefield definition
#-------------------------------------------------------------------------------------------------------------

# Define the probe LJ parameters (H2 for the time being--change here)
eps1  = 11.25 # in kelvin
sig1  = 2.68  # in angstrom
rcut = 14     # Maximum range under consideration for the LJ interactions


# Define grid points in a unit box (Non adaptive grid)
# Number of grid points
nx=20
ny=20
nz=20

# Intially the grids are defined on a unit box
x_grid = np.linspace(0,1,nx)
y_grid = np.linspace(0,1,ny)
z_grid = np.linspace(0,1,nz)

# Read the forcefield information from RASPA force field directory
# I made some changes there to make life easy
force_field = 'GenericMOFs'
#pa_file = open('./forcefield/'+forcefield+'/pseudo_atoms.def').readlines()
#NumberOfPseudoAtoms = int(pa_file.[1].split()[0])

ff_file = open('./forcefield/'+force_field+'/force_field_mixing_rules.def').readlines()

# --------------------------------------------------------------------------------------------------------------------
# Prepare dir structure
# ---------------------------------------------
path_orig = os.getcwd()
os.chdir('CIF_FILES')
path_files=os.getcwd()
cif_list=os.listdir('.')

#-------------------------------------------------------------------------------------------------------------
# The heart of the code
#-------------------------------------------------------------------------------------------------------------
for name_index in range(len(cif_list)):
	os.chdir(cif_list[name_index])
	cif_file_name = cif_list[name_index] + '.cif'
	
	
	# Read the coordinates from the cif file using pymatgen
	f1=CifParser(cif_file_name)
	struct=f1.get_structures()[0]
	
	
	# Need the box lengths to ensure the distance criterion
	# Create the cell matrix to ensure this criterion
	
	lx = struct.lattice.a
	ly = struct.lattice.b
	lz = struct.lattice.c
	
	alpha = struct.lattice.alpha * (np.pi/180.0)
	beta  = struct.lattice.beta * (np.pi/180.0)
	gamma = struct.lattice.gamma * (np.pi/180.0)
	
	
	
	unit_cell_matrix = [[1 , 0, 0],[np.cos(gamma),np.sin(gamma),0],[0, 0, 1]]
	unit_rec_space = np.linalg.inv(unit_cell_matrix)
	
	
	# Make super cell based on the cut off criteria
	# the minimum box edge must be greater than twice the maximum range of interaction  
	# set to fixed 14 angstroms
	
	nx_cells = np.ceil(2*rcut/lx)    # magic formula
	ny_cells = np.ceil(2*rcut/ly)
	nz_cells = np.ceil(2*rcut/lz) 
	
	
	struct.make_supercell([nx_cells,ny_cells,nz_cells])
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
	
	
	# Read the corresponding forcefield parameters (pseudo atoms information)
	# for the atoms in the MOF
	
	#Get atom names in the MOF
	mof_atm_names = []
	for i in range(Number_Of_Atoms):
		mof_atm_names.append(str(struct.species[i]))
	
	
	# Get the atom indices from the force_field file for the atoms in the MOF
	mof_atm_indices = np.zeros((Number_Of_Atoms,1))
	for i in range(len(mof_atm_names)):
		for line in ff_file:
			if mof_atm_names[i] == line.split()[0]:
				mof_atm_indices[i]=ff_file.index(line) 
	
	pot=np.zeros((nx,ny,nz))
	
	#Calculate LJ interaction energy
	for atm_index in range(Number_Of_Atoms):
		# Lorentz Berthelot
		sig2 = float(ff_file[int(mof_atm_indices[atm_index][0])].split()[3])  # Get sigma
		eps2 = float(ff_file[int(mof_atm_indices[atm_index][0])].split()[2])  # get epsilon

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

					#Transform back to find the actual distance
					r_act = np.dot(A , drij)
					rsq=r_act[0]**2+r_act[1]**2+r_act[2]**2

					if np.sqrt(rsq) <= rcut:
						pot[i][j][k] += lj(eps,sig,rsq)
	
	#-------------------------------------------------------------------------------------------------------------
	# Output the VTS, grid energy values and attractive zone
	#-------------------------------------------------------------------------------------------------------------

	N_grid_total = nx *ny *nz
	
	#Write the VTK file
	# Define the crazy unstructured grid      
	dx, dy, dz = 1.0/nx, 1.0/ny, 1.0/nz     
	X = np.arange(0, 1 + 0.1*dx, dx, dtype='float64') 
	Y = np.arange(0, 1 + 0.1*dy, dy, dtype='float64') 
	Z = np.arange(0, 1 + 0.1*dz, dz, dtype='float64') 
	
	x = np.zeros((nx , ny , nz)) 
	y = np.zeros((nx , ny , nz)) 
	z = np.zeros((nx , ny , nz)) 
	
	for k in range(nz): 
		for j in range(ny):
			for i in range(nx): 
				x[i,j,k] = X[i] 
				y[i,j,k] = Y[j] 
				z[i,j,k] = Z[k]
	
	for k in range(nz):
		for j in range(ny):
			for i in range(nx):
				[x[i,j,k], y[i,j,k],z[i,j,k]] = np.dot(A,[x[i,j,k], y[i,j,k],z[i,j,k]])
	
	# Write pot into .vts file
	gridToVTK("./pot", x, y, z, pointData = {"Potential" : pot})             
	
	
	# Histogram into bins (predefined?)
	e_vals = np.reshape(pot,(N_grid_total,1))  # Reshape into linear array
	bins1 = np.linspace(min(e_vals),0,31)
	e_hist, binedges1 = np.histogram(e_vals,bins=bins1, normed = 'false') # Histogram
	bincenters = 0.5*(binedges1[1:]+binedges1[:-1]) # Bincenters
	data = np.vstack((bincenters, e_hist) )
	np.savetxt('histogram.txt',data) # Write histogram data
	 
	 
	#Write the raw energy values        
	f3=open('Details.txt','w')
	f3.write(str(nx)+'\t'+str(ny)+'\t'+str(nz)+'\n')
	f3.write(str(nx_cells)+'\t'+str(ny_cells)+'\t'+str(nz_cells)+'\n')
	f3.write(str(lx)+'\t'+str(ly)+'\t'+str(lz)+'\n')       
	f3.write(str(alpha)+'\t'+str(beta)+'\t'+str(gamma)+'\n')       
	f3.close()
	
	f3=open('Energy_Values.txt','a')                
	np.savetxt(f3,e_vals)
	f3.close()
	
	
	#Print the fraction of the attractive zone
	f2=open('../../Metric.txt','a')
	f2.write('The attraction zone for '+cif_file_name+':\t')
	f2.write(str(float(sum((e_vals < 0) & (e_vals>min(e_vals)))[0])/N_grid_total) + '\n')
	f2.close()
	
	# write the raw grid data into text file straight array
	# Write the corresponding xyz coordinates
	# Write the corresponding MOF coordinates to .xyz file using pymatgen
	f4=open('Structure_Information.dat','w')
	f4.write(str(struct))
	f4.close()
	
	os.chdir(path_files)
os.chdir(path_orig)