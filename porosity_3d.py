from __future__ import division
import numpy as np
import math

'physics'
Ly = 25      #domain length in y, m
Lx = 10      #domain length in x, m
Lz = 10       #domain length in z, m
k0 = 2       #permeability, m2
mu = 1       #viscosity of the fluid, kg/(ms)
eta = 1       # effective viscosity
B = 0.1     #compressibility
g = -1       #acceleration due to gravity, m/s2
drho = 1       #density of the fluid, kg/m3
phi_0 = 0       #porosity

'numerics'      #play with the resulation to see if we still have the same results
ny = 128         # number of grid points in y
nx = int((ny*Lx/Ly))  # number of grid points in x
nz = int((ny*Lz/Ly))  # number of grid points in z
dt = 0.00005        # time step 
nt = 1000000        # number of time steps 
nk = 10000      #number of sub-iterations

'preprocessing'
dx = Lx/(nx-1)                # grid spacing in x
dy = Ly/(ny-1)                # grid spacing in y
dz = Ly/(nz-1)                # grid spacing in z
x = np.zeros([nx,ny,nz])
y = np.zeros([nx,ny,nz])
z = np.zeros([nx,ny,nz])

for ix in xrange (0,nx):
    for iy in xrange (0,ny):
        for iz in xrange (0,nz):
            x[ix,iy,iz] = (ix)*dx
            y[ix,iy,iz] = (iy)*dy
            z[ix,iy,iz] = (iz)*dz

K_m = k0/mu

'initial'
P = np.zeros(x.shape)           #pressure
P_n = np.zeros(x.shape)
P_old = np.zeros(x.shape)
res = np.zeros(x.shape)
phi = phi_0 + np.exp((-( (x-Lx*0.5)/2.5 )**2 -( (y-Ly*0.1)/2.5)**2 -( (z - Lz*0.5)/2.5 )**2 ))
# playing with 0.1 and 0.5 will change the appearance

for j in xrange(int(5/dy)+1,int(7/dy)+1):
    for i in xrange(0,nx):
        for k in xrange(0,nz):
            phi[i,j,k] = 0.1

drhog = drho*g
K = np.multiply(K_m, (np.power (phi,3)))
K_avx = np.multiply((K[:-1,:,:] + K[1:,:,:]),0.5)
K_avy = np.multiply((K[:,:-1,:] + K[:,1:,:]),0.5)
K_avz = np.multiply((K[:,:,:-1] + K[:,:,1:]),0.5)
P[0,:,:] = 0
P[-1,:,:] = 0
P[:,0,:] = 0
P[:,-1,:] = 0
P[:,:,0] = 0
P[:,:,-1] = 0
qx = np.multiply((-K_avx), ((P[1:,:,:]-P[:-1,:,:])/dx))
qy = np.multiply((-K_avy), (((P[:,1:,:]-P[:,:-1,:])/dy)+drhog))
qz = np.multiply((-K_avz), ((P[:,:,1:]-P[:,:,:-1])/dz))
dPdt = np.multiply((-1/B),((qx[1:,1:-1,1:-1]-qx[:-1,1:-1,1:-1])/dx + (qy[1:-1,1:,1:-1]-qy[1:-1,:-1,1:-1])/dy + (qz[1:-1,1:-1,1:]-qz[1:-1,1:-1,:-1])/dz + P[1:-1,1:-1,1:-1]/eta))
P[1:-1,1:-1,1:-1] = P[1:-1,1:-1,1:-1] + dt*dPdt


'action'
for it in xrange(0,nt):
        print it
        P_n[:,:,:] = P[:,:,:]
        for ik in xrange(0,nk):
                print ik
                P_old[:,:,:] = P[:,:,:]
                qx = np.multiply((-K_avx),((P[1:,:,:]-P[:-1,:,:])/dx))
                qy = np.multiply((-K_avy), (((P[:,1:,:]-P[:,:-1,:])/dy)+drhog))
                qz = np.multiply((-K_avz),((P[:,:,1:]-P[:,:,:-1])/dz))
                dPdt = np.multiply((-1/B),((qx[1:,1:-1,1:-1]-qx[:-1,1:-1,1:-1])/dx + (qy[1:-1,1:,1:-1]-qy[1:-1,:-1,1:-1])/dy + (qz[1:-1,1:-1,1:]-qz[1:-1,1:-1,:-1])/dz + P[1:-1,1:-1,1:-1]/eta))
                P[1:-1,1:-1,1:-1] = P_n[1:-1,1:-1,1:-1] + dt*dPdt
                res = np.fabs(np.fabs(P_old) - np.fabs(P))
                res_max = res.max()
                res_total = res.sum(0).sum(0).sum(0)
                P_old_total = (np.fabs(P_old)).sum(0).sum(0).sum(0)
                res_avg = np.divide(res_total, P_old_total)

                print res_max
                if res_max < 0.000001:
                        break

        phi[1:-1,1:-1,1:-1] = phi[1:-1,1:-1,1:-1] + dt*(B*dPdt + P[1:-1,1:-1,1:-1]/eta)
