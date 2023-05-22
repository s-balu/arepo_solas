import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import matplotlib.pyplot as plt
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d
from scipy.optimize import root

def rho(x):
    M = 10**10
    a = 7670
    fg = 0.16
    
    return fg*M/(2*np.pi) * a/x * 1/(x+a)**3


def cdf(x):
    integral, _ = quad(rho, 0.1, x)
    normalization, _ = quad(rho, 0.1, 4)
    return integral / normalization

def sample_from_pdf(num_samples):
    samples = []
    for _ in range(num_samples):
        u = np.random.uniform(0, 1)
        # Find the x value that corresponds to the given u
        x = root(lambda x: cdf(x) - u, x0=0.1, method='hybr').x[0]
        samples.append(x)
        print(len(samples))
    return samples

def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z
    
        
def dphi_dr(x):
    G = 4.301871e-03
    M = 10**10
    a = 7670

    if(x > 0):
      return G*M/(x+a)**2
    else:
      print("r is zero")
      return 0
      
def hydro_eqn(p, x, rho, dphi_dr):
 
    return -rho(x) * dphi_dr(x)


#def dphi_dr(rt):
#    if(rt > 0):
#      return np.log(1+rt)/(rt**2) - 1/(1+rt)/rt
#    else:
#      print("r is zero")
#      return 0
   
simulation_directory = str(sys.argv[1])
print("BHs: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

# double precision: np.float64, for single use np.float32
FloatType = np.float64  
IntType = np.int32

#blackhole ics
PosBh = np.zeros(3, dtype=FloatType)
PosBh[0] = 4
PosBh[1] = 4
PosBh[2] = 4

MassBh = 10**6
hsml   = 0.05

#generate random samples from the distribution (i.e. HQ rho)
N=10**4
R = sample_from_pdf(N)
print("Sample done")
#convert the radial distances to Cartesian coordinates
x_values = []
y_values = []
z_values = []

for r in R:
    #generate random values for theta and phi
    phi = np.arccos(np.random.uniform(-1, 1))    # Cosine of Angle between -1 and 1 (polar angle)
    theta = np.random.uniform(0, 2*np.pi)    # Angle between 0 and 2*pi (azimuthal angle)

    #convert from spherical to Cartesian coordinates
    x, y, z = spherical_to_cartesian(r,phi,theta)

    #append the coordinates to the respective lists
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)

x_values = np.array(x_values, dtype=FloatType)
y_values = np.array(y_values, dtype=FloatType)
z_values = np.array(z_values, dtype=FloatType)
print("Conversion done")
Boxsize = FloatType(8)
"""CellsPerDimension = IntType(100)
NumberOfCells = CellsPerDimension * CellsPerDimension * CellsPerDimension

dx = Boxsize / FloatType(CellsPerDimension)
pos_first, pos_last = FloatType(0.5) * dx, Boxsize - FloatType(0.5) * dx

Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = xx.reshape(NumberOfCells)
Pos[:,1] = yy.reshape(NumberOfCells)
Pos[:,2] = zz.reshape(NumberOfCells)
for i in range(NumberOfCells):
    Pos[i,0]+=np.random.uniform(low=-0.5*dx,high=0.5*dx)
    Pos[i,1]+=np.random.uniform(low=-0.5*dx,high=0.5*dx)
    Pos[i,2]+=np.random.uniform(low=-0.5*dx,high=0.5*dx)

#sort by distance from bh
x = Pos[:,0] - PosBh[0]
y = Pos[:,1] - PosBh[1]
z = Pos[:,2] - PosBh[2]
r = np.sqrt(x**2 + y**2 + z**2)

for i in range(N):
    if(R[i] == 0):
        r = np.delete(r, i)
        x_values = np.delete(x, i)
        y_values = np.delete(y, i)
        z = np.delete(z, i)
"""        

r_sort    = np.sort(R)
r_reverse = r_sort[::-1]

#same mass for sph particles
Mass = np.full(N, 4.35*10**(-4), dtype=FloatType)

#find pressure from hydrostatic equilibrium
p0 = 0
sol_reverse = odeint(hydro_eqn, p0, r_reverse, args=(rho,dphi_dr))
sol         = sol_reverse[::-1]
Pressure    = interp1d(r_sort, sol[:,0], kind='cubic')(r_sort)

gamma = FloatType(5/3) 
gamma_minus1 = gamma - FloatType(1.0)

#utherm for perfect gas
Density=np.zeros(N,dtype=FloatType)
for i in range(N):
    Density[i] = rho(r_sort[i])

Uthermal = Pressure / gamma_minus1 / Density

#sort positions
indices = np.argsort(R)
x_sorted = x_values[indices]
y_sorted = y_values[indices]
z_sorted = z_values[indices]
Pos_sorted = np.zeros([N,3], dtype=FloatType)

Pos_sorted[:,0] = x_sorted
Pos_sorted[:,1] = y_sorted
Pos_sorted[:,2] = z_sorted

Pos_sorted[:,0] += PosBh[0]
Pos_sorted[:,1] += PosBh[1]
Pos_sorted[:,2] += PosBh[2]

print(min(Pos_sorted[:,2]), max(Pos_sorted[:,2]))

bins = 100
count, bin_edges = np.histogram(r_sort, bins=bins)
bin_edges = np.delete(bin_edges, -1)
count = np.array(count, dtype=FloatType)
count/=N

plt.scatter(bin_edges, count)
plt.plot(bin_edges, rho(bin_edges)/435)
plt.show()


""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells
part5 = IC.create_group("PartType5")  

## write header entries
NumPart = np.array([N, 0, 0, 0, 0, 1], dtype=IntType)
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
if Pos_sorted.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, N+1))
part0.create_dataset("Coordinates", data=Pos_sorted)
part0.create_dataset("Masses", data=Mass)
#part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)
part5.create_dataset("ParticleIDs" , data = np.arange(N+1,N+2))
part5.create_dataset("Coordinates", data=PosBh)
part5.create_dataset("Masses" , data = MassBh)
#part5.create_dataset("Velocities", data=VelBh[0:1])
part5.create_dataset("BlackholeHsml" , data = hsml)



## close file
IC.close()

""" normal exit """
sys.exit(0) 
