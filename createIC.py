import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
from scipy.integrate import odeint
from scipy.interpolate import interp1d


#def dphi_dr(rt):
#    if(rt > 0):
#      return np.log(1+rt)/(rt**2) - 1/(1+rt)/rt
#    else:
#      print("r is zero")
#      return 0
   
   
def rho(rt):
    M = 10**10
    a = 7670
    fg = 0.16
    
    return fg*M/(2*np.pi) * a/rt * 1/(rt+a)**3
    
       
def dphi_dr(rt):
    G = 4.301871e-03
    M = 10**10
    a = 7670

    if(rt > 0):
      return G*M/(rt+a)**2
    else:
      print("r is zero")
      return 0
      
      
def hydro_eqn(p, rt, rho, dphi_dr):
 
    return -rho(rt) * dphi_dr(rt)


simulation_directory = str(sys.argv[1])
print("BHs: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

#blackhole ics
PosBh = np.zeros(3, dtype=FloatType)
PosBh[0] = 2.5
PosBh[1] = 2.5
PosBh[2] = 2.5

MassBh = 10**6
hsml   = 0.05

#make grid
Boxsize = FloatType(5)
CellsPerDimension = IntType(100)
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

for i in range(NumberOfCells):
    if(r[i] == 0):
        r = np.delete(r, i)
        x = np.delete(x, i)
        y = np.delete(y, i)
        z = np.delete(z, i)
        
indices = np.argsort(r)

x_sorted = x[indices]
y_sorted = y[indices]
z_sorted = z[indices]

Pos_sorted = np.zeros([NumberOfCells,3], dtype=FloatType)
Pos_sorted[:,0] = x_sorted
Pos_sorted[:,1] = y_sorted
Pos_sorted[:,2] = z_sorted

r_sort    = np.sort(r)
r_sorted  = np.unique(r)
r_reverse = r_sorted[::-1]

#get density from hernquist formula
Density = np.zeros(NumberOfCells, dtype=FloatType)
for i in range(NumberOfCells):
    Density[i] = rho(r_sort[i])
    
Volume      = np.full(NumberOfCells, dx*dx*dx, dtype=FloatType)
Mass        = Density * Volume
print(np.mean(Mass)) #simulating r < 4pc

#find pressure from hydrostatic equilibrium
p0 = 1

sol_reverse = odeint(hydro_eqn, p0, r_reverse, args=(rho,dphi_dr))
sol         = sol_reverse[::-1]
Pressure    = interp1d(r_sorted, sol[:,0], kind='cubic')(r_sort)

gamma = FloatType(5/3)  ## note: this has to be consistent with the parameter settings for Arepo!
gamma_minus1 = gamma - FloatType(1.0)

Uthermal = Pressure / gamma_minus1 / Density

Pos_sorted[:,0] += PosBh[0]
Pos_sorted[:,1] += PosBh[1]
Pos_sorted[:,2] += PosBh[2]

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells
part5 = IC.create_group("PartType5")  

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 1], dtype=IntType)
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

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1))
part0.create_dataset("Coordinates", data=Pos_sorted)
part0.create_dataset("Masses", data=Mass)
#part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)
part5.create_dataset("ParticleIDs" , data = np.arange(NumberOfCells+1,NumberOfCells+2))
part5.create_dataset("Coordinates", data=PosBh)
part5.create_dataset("Masses" , data = MassBh)
#part5.create_dataset("Velocities", data=VelBh[0:1])
part5.create_dataset("BlackholeHsml" , data = hsml)



## close file
IC.close()

""" normal exit """
sys.exit(0) 
