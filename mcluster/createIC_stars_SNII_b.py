import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import matplotlib.pyplot as plt
from scipy.integrate import  quad
from scipy.optimize import root

UnitMass = 1.989e33
UnitLength = 3.086e18
UnitVel = 1.e4
UnitTime = UnitLength/UnitVel
UnitDensity = UnitMass/UnitLength**3.
UnitEnergy = UnitMass*UnitVel**2.
UnitPressure = UnitEnergy / UnitLength**3.

Grav = 6.6742e-8
CLIGHT = 3.e10
PROTONMASS = 1.6726e-24
BOLTZMANN = 1.38065e-16
THOMPSON = 6.65e-25
PI = np.pi
HYDROGEN_MASSFRAC = 0.6
GAMMA_MINUS1 = 5./3. - 1.
PC = 3.086e18

def rho(x):
    M = 10**10
    a = 7670
    fg = 0.16
    
    return 4*np.pi*x**2 * fg*M/(2*np.pi)*a/x*1/(x+a)**3 #0.2

def rho_s(x):
    M  = 2.24*10**6
    a  = 4
    fg = 1
    
    return 4*np.pi*x**2 * fg*M/(2*np.pi)*a/x*1/(x+a)**3
    
def rho_b(x):
    M = 10**5
    a = 7670
    fg = 0.16
    
    return 4*np.pi*x**2 * fg*M/(2*np.pi)*a/x*1/(x+a)**3 #2*10**(-6)

def cdf(x, flag):
    if(flag == 0):
        integral, _ = quad(rho, 0, x)
        normalization, _ = quad(rho, 0, 40)
    if(flag == 1):
        integral, _ = quad(rho_s, 0, x)
        normalization, _ = quad(rho_s, 0, 40)
    if(flag == 2):
        integral, _ = quad(rho_b, 40, x)
        normalization, _ = quad(rho_b, 40, 80)   
    
    return integral / normalization

def sample_from_pdf(num_samples, flag):
    samples = []
    for _ in range(num_samples):
        u = np.random.uniform(0, 1)
        # Find the x value that corresponds to the given u
        x = root(lambda x: cdf(x,flag) - u, x0=0.1, method='hybr').x[0]
        samples.append(x)
        print(len(samples))
    return samples

def spherical_to_cartesian(r,phi,theta):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

def density(x):
    M  = 10**10
    a  = 7670
    fg = 0.16
    
    return fg*M/(2*np.pi)*a/x*1/(x+a)**3

def density_s(x):
    M  = 2.24*10**6
    a  = 4
    fg = 1
    
    return fg*M/(2*np.pi)*a/x*1/(x+a)**3
    
def density_b(x):
    M  = 10**5
    a  = 7670
    fg = 0.16
    
    return fg*M/(2*np.pi)*a/x*1/(x+a)**3

def dm(x, flag):
    if(flag == 0):
        return 4*np.pi*x**2 * density(x)
    if(flag == 1):
        return 4*np.pi*x**2 * density_s(x)
    if(flag == 2):
        return 4*np.pi*x**2 * density_b(x)

#Generate ICs
simulation_directory = str(sys.argv[1])
print("BHs: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

#double precision: np.float64, for single use np.float32
FloatType = np.float64  
IntType = np.int32

#Blackhole ICs
PosBh = np.zeros(3, dtype=FloatType)
PosBh[0] = 80
PosBh[1] = 80
PosBh[2] = 80

hsml   = 5

#Generate random samples from the distribution (i.e. HQ rho)
N = 10**6
Nb  = 10**5
R = sample_from_pdf(N,0)
Rb = sample_from_pdf(Nb,2)

R  = np.array(R)
Rb = np.array(Rb)
print(R.min(), R.max(), Rb.min(), Rb.max())
Ntot = N + Nb
R = np.append(R,Rb)

print("Sample done")
#convert the radial distances to Cartesian coordinates
x_values = []
y_values = []
z_values = []

for r in R:
    #generate random values for theta and phi
    phi = np.arccos(1-2*np.random.uniform(0, 1)) # Cosine of Angle between -1 and 1 (polar angle)
    theta = 2*np.pi * np.random.uniform(0, 1)    # Angle between 0 and 2*pi (azimuthal angle)

    #convert from spherical to Cartesian coordinates
    x, y, z = spherical_to_cartesian(r,phi,theta)

    #append the coordinates to the respective lists
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)

R        = np.array(R, dtype=FloatType)
x_values = np.array(x_values, dtype=FloatType)
y_values = np.array(y_values, dtype=FloatType)
z_values = np.array(z_values, dtype=FloatType)

print("Conversion done")

#sort gas cell arrays
indices = np.argsort(R)        
r = R[indices]
x = x_values[indices]
y = y_values[indices]
z = z_values[indices]

#insert positions
Boxsize = FloatType(160)
Pos = np.zeros([Ntot,3],dtype=FloatType)

Pos[:,0] = x
Pos[:,1] = y
Pos[:,2] = z

Pos[:,0] += PosBh[0]
Pos[:,1] += PosBh[1]
Pos[:,2] += PosBh[2]

#same mass for sph particles
M, _ = quad(dm, 0, 40, args=(0))
m = M/N
Mass = np.full(N, m, dtype=FloatType)
print("m:", m)

Mb, _ = quad(dm, 40, 80, args=(2))
mb = Mb/Nb
Massb = np.full(Nb, mb, dtype=FloatType)
print("mb:", mb)

Mass = np.append(Mass, Massb)

print("Mass:", Mass)

#pressure
Pressure = np.zeros(Ntot, dtype=FloatType)

T = 10**3
utherm_0  = T * BOLTZMANN / GAMMA_MINUS1 / PROTONMASS / 0.6
utherm_0 /= (UnitEnergy/UnitMass)

for i in range(Ntot):
    if r[i] > 40:
        Pressure[i] = GAMMA_MINUS1*utherm_0*density(r[i])
    else:
        Pressure[i] = GAMMA_MINUS1*utherm_0*density(r[i])
    
#utherm
Density=np.zeros(Ntot,dtype=FloatType)

for i in range(Ntot):
    if r[i] > 40:
        Density[i] = density_b(r[i])
    else:
        Density[i] = density(r[i])

Uthermal = Pressure / GAMMA_MINUS1 / Density

d = Density[r<40]
mean = np.mean(d)

d_b = Density[r>40]
mean_b = np.mean(d_b)

print("AVG Densities:", mean, mean_b)

#checks
bins = 100
count, bin_edges = np.histogram(r, bins=bins)
bin_edges = np.delete(bin_edges, -1)
bin_width = bin_edges[1]-bin_edges[0]
count = np.array(count, dtype=FloatType)
count/= (N * bin_width)

rho_analytic = np.zeros(bins)
for i in range(bins):
    rho_analytic[i] = rho(bin_edges[i])
norm, _ = quad(rho, 0, 40)
rho_analytic /= norm

plt.scatter(bin_edges, count)
plt.plot(bin_edges,rho_analytic)
plt.xlim(0,40)
plt.show()

bins = 100
count, bin_edges = np.histogram(r, bins=bins)
bin_edges = np.delete(bin_edges, -1)
bin_width = bin_edges[1]-bin_edges[0]
count = np.array(count, dtype=FloatType)
count/= (Nb * bin_width)

rho_analytic = np.zeros(bins)
for i in range(bins):
    rho_analytic[i] = rho_b(bin_edges[i])
norm, _ = quad(rho_b, 40, 80)
rho_analytic /= norm

plt.scatter(bin_edges, count)
plt.plot(bin_edges,rho_analytic)
plt.xlim(40,80)
plt.show()

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells
part5 = IC.create_group("PartType5")
#part4 = IC.create_group("PartType4")  

## write header entries
NumPart = np.array([Ntot,0,0,0,0,1], dtype=IntType)
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
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

Pos_star5  = [80,80,80]
Mass_star5 = 20
Vel_star5  = [0,0,0]   

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, Ntot+1))
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
#part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)
part5.create_dataset("ParticleIDs" , data=np.arange(Ntot+1, Ntot+2))
part5.create_dataset("Coordinates", data=Pos_star5)
part5.create_dataset("Masses" , data=Mass_star5)
part5.create_dataset("Velocities", data=Vel_star5)
part5.create_dataset("BlackholeHsml" , data=np.full(1, hsml, dtype=FloatType))
#part4.create_dataset("ParticleIDs" , data=np.arange(Ntot+1,Ntot+N4+1))
#part4.create_dataset("Coordinates", data=Pos_star4)
#part4.create_dataset("Masses" , data=Mass_star4)
#part4.create_dataset("Velocities", data=Vel_star4)

## close file
IC.close()

""" normal exit """
sys.exit(0) 
