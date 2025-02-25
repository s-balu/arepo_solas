import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import matplotlib.pyplot as plt
from scipy.integrate import  quad
from scipy.optimize import root

num_rows = 0

#Load starcluster
data = []
with open('test.txt', 'r') as file:
    for line in file:
        num_rows += 1
        values = line.strip().split()
        data.append([float(value) for value in values])
        

#separate the columns
ms = [row[0] for row in data]
xs = [row[1] for row in data]
ys = [row[2] for row in data]
zs = [row[3] for row in data]
vxs = [row[4] for row in data]
vys = [row[5] for row in data]
vzs = [row[6] for row in data]

#Generate ICs
simulation_directory = str(sys.argv[1])
print("BHs: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

#double precision: np.float64, for single use np.float32
FloatType = np.float64  
IntType = np.int32

Nstar = num_rows

Boxsize = FloatType(160)
Pos_star = np.zeros([Nstar,3],dtype=FloatType)

Pos_star[:,0] = xs
Pos_star[:,1] = ys
Pos_star[:,2] = zs

#Pos_star[:,0] += PosBh[0]
#Pos_star[:,1] += PosBh[1]
#Pos_star[:,2] += PosBh[2]

Mass_star = np.array(ms)
for i in range(Nstar):
    if ms[i] > 5:
        print(ms[i])

Vel_star = np.zeros([Nstar,3],dtype=FloatType)
Vel_star[:,0] = vxs
Vel_star[:,1] = vys
Vel_star[:,2] = vzs

from mpl_toolkits.mplot3d import Axes3D

#condition = Mass_star > 80

#xs = np.array(xs)[condition]
#ys = np.array(ys)[condition]
#zs = np.array(zs)[condition]
# Plot the particles on a sphere
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs, ys, zs)
#plt.show()

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part5 = IC.create_group("PartType5")  

## write header entries
NumPart = np.array([0,0,0,0,0,Nstar], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos_star.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part5.create_dataset("ParticleIDs" , data = np.arange(1, Nstar+1))
part5.create_dataset("Coordinates", data=Pos_star)
part5.create_dataset("Masses" , data = Mass_star)
part5.create_dataset("Velocities", data=Vel_star)
#part5.create_dataset("BlackholeHsml" , data = np.full(Nstar, hsml, dtype=FloatType))



## close file
IC.close()

""" normal exit """
sys.exit(0) 
