from __future__ import division
import numpy as np
import re
from scipy import interpolate
from scipy.interpolate import griddata

# Constant

hbar = 6.62607004e-34 # J S
pi = np.pi
Ry2J = 2.1798741e-18 # J/Ry
Bohr2m = 5.29177249e-11 # m/Bohr
me = 9.10938356e-31 # kg
latt = 10.458182 * Bohr2m


# The first Brillouin zone of a face centered cubic lattice

b1 = [  -1  ,  -1,     1]
b2 = [ 1   ,  1 ,    1]
b3 = [-1   ,  1 ,   -1]

# Transforming matrix (Crystal to Cartesian)

bm = np.matrix([b1,b2,b3])

# The band minimum

km = [0.0, 0.0, 0.0]

# Generating mesh

## Mesh dimensions

Nx = 5 # 5 is chosen for convenience of decrete derivative
Ny = 5
Nz = 5

pos = np.zeros((Nx*Ny*Nz,3))

## Resolution

dx = 1.0/100.0/(Nx-1)
dy = 1.0/100.0/(Ny-1)
dz = 1.0/100.0/(Nz-1)

# K points

for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            pos[i*Ny*Nz + j*Ny + k, 0] = km[0] + i * dx - dx * (Nx-1)/2
            pos[i*Ny*Nz + j*Ny + k, 1] = km[1] + j * dy - dy * (Ny-1)/2
            pos[i*Ny*Nz + j*Ny + k, 2] = km[2] + k * dz - dz * (Nz-1)/2

# Write the K points

f = open("kcdt.dat", "wb")
f.write(str(Nx * Ny * Nz)+'\n')
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            f.write(('{0:12.8f} {1:12.8f} {2:12.8f} {3:12.6e}'.format(pos[i*Ny*Nz + j*Ny + k, 0],
            pos[i*Ny*Nz + j*Ny + k, 1],
            pos[i*Ny*Nz + j*Ny + k, 2],
            1/(Nx * Ny * Nz)))+'\n')
f.close()

# Read the band structure



f = open("fit_wo_SO.out", "r") # output file for band calculation
lines = f.readlines()
for i, line in enumerate(lines):
    if re.search("End of band structure calculation", line):
       mark_line_start = i # find the start
    if re.search("Writing output data file", line):
       mark_line_end = i # find the end

# Number of bands

nbnd = 12
nline = int(round(nbnd/8.0)) # number of lines associated with k coordinates
spacing = 2*nline+5

kmat = np.zeros((Nx*Ny*Nz,3))
kmat_dxyz = np.zeros((Nx*Ny*Nz,3))
Emat = np.zeros((Nx*Ny*Nz,1))
Emat_dE = np.zeros((Nx*Ny*Nz,1))
Eband = np.zeros((Nx*Ny*Nz,nbnd))


# Define the index of band minimum (from low energy to high energy)

Nmin = 3 # 10 means No.10 band with the count starting from 1

# Read the k coordinates and energy

for i in range(Nx*Ny*Nz):
    for nnn, a in enumerate(lines[mark_line_start + spacing * i + 2].split()):
        if "=" in a:
            # lines[mark_line_start + spacing * i + 2].split()[nnn] = temp[1:]
            kmat[i,:] = np.dot(pos[i,:],bm) # float(temp[1:])
            # kmat[i,0] = float(lines[mark_line_start + spacing * i + 2].split()[1])
            # kmat[i,1] = # float(lines[mark_line_start + spacing * i + 2].split()[2])
            # kmat[i,2] = # float(lines[mark_line_start + spacing * i + 2].split()[3])
            for j in range(nline):
                for k in range(8):
                    if j*8+k < nbnd:
                        Eband[i,j*8+k] = float(lines[mark_line_start + spacing * i + j +4].split()[k])
    Emat[i,0] = Eband[i,Nmin -1]
f.close()

# Calculate the deviation

for i in range(Nx*Ny*Nz):
    kmat_dxyz[i,:] = kmat[i,:]-np.dot(km,bm)
    Emat_dE[i,:] = (Emat[i,:]-min(Emat))

# Interplate

temp1 =  np.dot([dx, dy ,dz],bm)


ddx = abs(temp1[0, 0])
ddy = abs(temp1[0, 1])
ddz = abs(temp1[0, 2])

# fit fxx

f_2 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 0, 0), method='linear')
f_1  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 0, 0), method='linear')
f0 = griddata(kmat_dxyz, Emat_dE, (0 * ddx, 0, 0), method='linear')
f1 = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 0, 0), method='linear')
f2 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 0, 0), method='linear')

fxx = 1.0/(12.0*ddx**2)*(-(f_2+f2)+16*(f_1+f1)-30*f0)

# fit fyy

f_2 = griddata(kmat_dxyz, Emat_dE, (0, -2 * ddy, 0), method='linear')
f_1  = griddata(kmat_dxyz, Emat_dE, (0, -1 * ddy, 0), method='linear')
f0 = griddata(kmat_dxyz, Emat_dE, (0, 0 * ddy, 0), method='linear')
f1 = griddata(kmat_dxyz, Emat_dE, (0, 1 * ddy, 0), method='linear')
f2 = griddata(kmat_dxyz, Emat_dE, (0, 2 * ddy, 0), method='linear')

fyy = 1.0/(12.0*ddy**2)*(-(f_2+f2)+16*(f_1+f1)-30*f0)


# fit fzz

f_2 = griddata(kmat_dxyz, Emat_dE, (0, 0, -2 * ddz), method='linear')
f_1  = griddata(kmat_dxyz, Emat_dE, (0, 0, -1 * ddz), method='linear')
f0 = griddata(kmat_dxyz, Emat_dE, (0, 0, 0 * ddz), method='linear')
f1 = griddata(kmat_dxyz, Emat_dE, (0, 0, 1 * ddz), method='linear')
f2 = griddata(kmat_dxyz, Emat_dE, (0, 0, 2 * ddz), method='linear')

fzz = 1.0/(12.0*ddz**2)*(-(f_2+f2)+16*(f_1+f1)-30*f0)

# fit fxy

f1_2 = griddata(kmat_dxyz, Emat_dE, (1 * ddx, -2 * ddy, 0), method='linear')
f2_1 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, -1 * ddy, 0), method='linear')
f_21 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 1 * ddy, 0), method='linear')
f_12 = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 2 * ddy, 0), method='linear')

f_1_2  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, -2 * ddy, 0), method='linear')
f_2_1  = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, -1 * ddy, 0), method='linear')
f12  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 2 * ddy, 0), method='linear')
f21  = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 1 * ddy, 0), method='linear')

f2_2 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, -2 * ddy, 0), method='linear')
f_22 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 2 * ddy, 0), method='linear')
f_2_2 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, -2 * ddy, 0), method='linear')
f22 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 2 * ddy, 0), method='linear')

f_1_1  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, -1 * ddy, 0), method='linear')
f11  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 1 * ddy, 0), method='linear')
f1_1  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, -1 * ddy, 0), method='linear')
f_11  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 1 * ddy, 0), method='linear')

fxy = 1.0/(600.0*ddx*ddy)*(-63*(f1_2+f2_1+f_21+f_12)\
    +63*(f_1_2+f_2_1+f12+f21)\
    +44*(f2_2+f_22-f_2_2-f22)\
    +74*(f_1_1+f11-f1_1-f_11))

# fit fxz

f1_2 = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 0, -2 * ddz), method='linear')
f2_1 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 0, -1 * ddz), method='linear')
f_21 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 0,  1 * ddz), method='linear')
f_12 = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 0,  2 * ddz), method='linear')

f_1_2  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 0,  -2 * ddz), method='linear')
f_2_1  = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 0,  -1 * ddz), method='linear')
f12  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 0,  2 * ddy), method='linear')
f21  = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 0,  1 * ddy), method='linear')

f2_2 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 0,  -2 * ddz), method='linear')
f_22 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 0,  2 * ddz), method='linear')
f_2_2 = griddata(kmat_dxyz, Emat_dE, (-2 * ddx, 0,  -2 * ddz), method='linear')
f22 = griddata(kmat_dxyz, Emat_dE, (2 * ddx, 0,  2 * ddz), method='linear')

f_1_1  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 0,  -1 * ddz), method='linear')
f11  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 0,  1 * ddz), method='linear')
f1_1  = griddata(kmat_dxyz, Emat_dE, (1 * ddx, 0,  -1 * ddz), method='linear')
f_11  = griddata(kmat_dxyz, Emat_dE, (-1 * ddx, 0,  1 * ddz), method='linear')

fxz = 1.0/(600.0*ddx*ddz)*(-63*(f1_2+f2_1+f_21+f_12)\
    +63*(f_1_2+f_2_1+f12+f21)\
    +44*(f2_2+f_22-f_2_2-f22)\
    +74*(f_1_1+f11-f1_1-f_11))

# fit fyz

f1_2 = griddata(kmat_dxyz, Emat_dE, (0, 1 * ddy,  -2 * ddz), method='linear')
f2_1 = griddata(kmat_dxyz, Emat_dE, (0, 2 * ddy,  -1 * ddz), method='linear')
f_21 = griddata(kmat_dxyz, Emat_dE, (0, -2 * ddy,   1 * ddz), method='linear')
f_12 = griddata(kmat_dxyz, Emat_dE, (0, -1 * ddy,  2 * ddz), method='linear')

f_1_2  = griddata(kmat_dxyz, Emat_dE, (0, -1 * ddy,  -2 * ddz), method='linear')
f_2_1  = griddata(kmat_dxyz, Emat_dE, (0, -2 * ddy,  -1 * ddz), method='linear')
f12  = griddata(kmat_dxyz, Emat_dE, (0, 1 * ddy,  2 * ddy), method='linear')
f21  = griddata(kmat_dxyz, Emat_dE, (0, 2 * ddy, 1 * ddy), method='linear')

f2_2 = griddata(kmat_dxyz, Emat_dE, (0, 2 * ddy,  -2 * ddz), method='linear')
f_22 = griddata(kmat_dxyz, Emat_dE, (0, -2 * ddy,  2 * ddz), method='linear')
f_2_2 = griddata(kmat_dxyz, Emat_dE, (0, -2 * ddy,  -2 * ddz), method='linear')
f22 = griddata(kmat_dxyz, Emat_dE, (0, 2 * ddy,  2 * ddz), method='linear')

f_1_1  = griddata(kmat_dxyz, Emat_dE, (0, -1 * ddy,  -1 * ddz), method='linear')
f11  = griddata(kmat_dxyz, Emat_dE, (0, 1 * ddy,  1 * ddz), method='linear')
f1_1  = griddata(kmat_dxyz, Emat_dE, (0, 1 * ddy,   -1 * ddz), method='linear')
f_11  = griddata(kmat_dxyz, Emat_dE, (0, -1 * ddy,   1 * ddz), method='linear')

fyz = 1.0/(600.0*ddy*ddz)*(-63*(f1_2+f2_1+f_21+f_12)\
    +63*(f_1_2+f_2_1+f12+f21)\
    +44*(f2_2+f_22-f_2_2-f22)\
    +74*(f_1_1+f11-f1_1-f_11))

# Unit conversion

cvt = Ry2J / (2.0 * pi/latt) ** 2 / hbar**2

# Inverse effective mass tensor

m_1 = cvt * \
np.matrix([[fxx[0],fxy[0],fxz[0]],[fxy[0],fyy[0],fyz[0]],[fxz[0],fyz[0],fzz[0]]])

# Effective mass tensor

m = np.linalg.inv(m_1)/me

# Eigenvalue and eigen vector (principle axies) of effective mass tensor

meig, eivtr = np.linalg.eigh(m, UPLO='L')

# Write the effective tensor

f = open("mass_tensor.out", "w") # output file for band calculation
f.write("Effective mass on priciple axies:\n\n")
f.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(meig[0],meig[1],meig[2])+'\n') # Eigenvector
f.write("\n")
f.write("Eigenvector:\n\n")
for i in range(3):
    f.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(eivtr[i,0],eivtr[i,1],eivtr[i,2])+'\n') # Eigenvector
f.write("\n")
f.write("Initial matrix:\n\n")
for i in range(3):
    f.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(m[i,0],m[i,1],m[i,2])+'\n') # Eigenvector
f.write("\n")
f.write("Conductivity effective mass (if cubic):\n\n")
mc = 1.0/3.0*(1/meig[0]+1/meig[1]+1/meig[2])
mc = 1/mc
f.write('{0:12.8f}'.format(mc)+'\n') # Eigenvector
md = (abs(meig[0])*abs(meig[1])*abs(meig[2]))**(1.0/3.0)
f.write("\n")
f.write("Density of state effective mass:\n\n")
f.write('{0:12.8f}'.format(md)+'\n') # Eigenvector
f.close()
