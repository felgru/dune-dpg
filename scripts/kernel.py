#!/usr/bin/python
# -*- coding: utf8 -*-
import numpy as np
from scipy.linalg import norm # to compute norm of an array
import matplotlib as mpl
import matplotlib.pyplot as plt

# Scaling function $\phi_{0,0}$ on the interval [-r,r[
def sf(x,r):
    normfactor=1./(np.sqrt(2*r))
    return np.array([normfactor*(val>=-r and val<r) for val in x])

# Wavelet $\psi_{n,k}$ on the interval [-r,r[
# Its support is in [xmin,xmax[ where
#     xmin = r*k*pow(2,-j+1)-r
#     xmax = r*(k+1)*pow(2,-j+1)-r
def wlt(j,k,x,r):
    xmin=r*k*pow(2,-j+1)-r
    xmax=r*(k+1)*pow(2,-j+1)-r
    xmiddle=(xmin+xmax)/2
    normfactor=1./(np.sqrt(xmax-xmin))
    positivePart=np.array([(val>=xmin and val<xmiddle) for val in x])
    negativePart=np.array([-1.*(val>=xmiddle and val<=xmax) for val in x])
    return normfactor*(positivePart+negativePart)

# Kernel
#     x in [-pi,pi]
#     y in [-pi,pi]
#     g in [0,1[
def phi(x,y,g):
    value=(1-g**2)/((2*np.pi)*(1+g**2-2*g*np.cos(x-y)))
    return value

def quadpoints(jx, kx, r):
    xmin=r*kx*pow(2,-jx+1)-r
    xmax=r*(kx+1)*pow(2,-jx+1)-r
    x=np.linspace(xmin,xmax,num=ngrid,endpoint=True)
    return (x, (xmax-xmin)/ngrid)

# ============
# Main program
# ============
ngrid=100 # Size mesh
J=7
r=np.pi
g=0.5 # forward peak coef

jy=0
ky=0

nwlt=0
result = []
(y, yquadweight) = quadpoints(0, 0, r)
# get wlt_{jy,ky}(y)
sfy=sf(y,r)
# print(np.sqrt((ymax-ymin)*np.mean(wlty*wlty))) # L2 norm wlt
# initialize resultx list
resultx = []
(x, xquadweight) = quadpoints(0, 0, r)
# get wlt_{jx,kx}(x)
sfx=sf(x,r)
# Values of the kernel
X, Y = np.meshgrid(x, y)
PHI=phi(X,Y,g)
# Product quadweight*kernel*wlt(x)*wlt(y)
quadweight = yquadweight * xquadweight
prod=quadweight*PHI*np.outer(sfx,sfy)
# Integral over [-pi,pi]x[-pi,pi]
resultx.append(prod.sum())
for jx in range(J):
    for kx in range(2**jx):
        (x, xquadweight) = quadpoints(jx, kx, r)
        # get wlt_{jx,kx}(x)
        wltx=wlt(jx,kx,x,r)
        # Values of the kernel
        X, Y = np.meshgrid(x, y)
        PHI=phi(X,Y,g)
        # Product quadweight*kernel*wlt(x)*wlt(y)
        quadweight = yquadweight * xquadweight
        prod=quadweight*PHI*np.outer(wltx,sfy)
        # Integral over [-pi,pi]x[-pi,pi]
        resultx.append(prod.sum())
result.append(resultx)
for jy in range(J):
    nwlt+=2**jy # This is just to count the number of wlts
    for ky in range(2**jy):
        (y, yquadweight) = quadpoints(jy, ky, r)
        # get wlt_{jy,ky}(y)
        wlty=wlt(jy,ky,y,r)
        # print(np.sqrt((ymax-ymin)*np.mean(wlty*wlty))) # L2 norm wlt
        # initialize resultx list
        resultx = []
        (x, xquadweight) = quadpoints(0, 0, r)
        # get wlt_{jx,kx}(x)
        sfx=sf(x,r)
        # Values of the kernel
        X, Y = np.meshgrid(x, y)
        PHI=phi(X,Y,g)
        # Product quadweight*kernel*wlt(x)*wlt(y)
        quadweight = yquadweight * xquadweight
        prod=quadweight*PHI*np.outer(sfx,wlty)
        # Integral over [-pi,pi]x[-pi,pi]
        resultx.append(prod.sum())
        for jx in range(J):
            for kx in range(2**jx):
                (x, xquadweight) = quadpoints(jx, kx, r)
                # get wlt_{jx,kx}(x)
                wltx=wlt(jx,kx,x,r)
                # Values of the kernel
                X, Y = np.meshgrid(x, y)
                PHI=phi(X,Y,g)
                # Product quadweight*kernel*wlt(x)*wlt(y)
                quadweight = yquadweight * xquadweight
                prod=quadweight*PHI*np.outer(wltx,wlty)
                # Integral over [-pi,pi]x[-pi,pi]
                resultx.append(prod.sum())
        result.append(resultx)

# ============
# Plot result
# ============
# # Plot with contourf
# xwlt=np.linspace(1,nwlt,nwlt)
# ywlt=xwlt
# Xwlt, Ywlt = np.meshgrid(xwlt, ywlt)
# origin = 'upper'
# CS = plt.contourf(Xwlt, Ywlt, result, 20,
#                   alpha=0.8,
#                   cmap=plt.cm.bone,
#                   origin=origin,
#                   interpolation='none')
# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = plt.colorbar(CS)
# # cbar.ax.set_ylabel('title of colorbar')


# Plot with mathshow
# mpl.rc('text',usetex=True)
MS = plt.matshow(result)
cbar = plt.colorbar(MS)
titleStr='Wlt analysis with $0\leq j\leq$'+str(J)+' for \n' \
    + r"$\phi(\theta,\theta')=(1-\gamma^2)/(1+\gamma^2-2\gamma\cos(\theta-\theta'))}$"+r'$,\ \gamma=$'+str(g)
plt.title(titleStr)
plt.xlabel('$(j,k)$')
plt.ylabel('$(j\',k\')$')
# plt.show()

titleSave="wlt_analysis-J_"+str(J)+"-gamma_"+str(g)+".pdf"
plt.savefig(titleSave, bbox_inches="tight")
plt.clf()

# plot logarithmic wavelet representation
log_result = []
for row in result:
    log_row = map(lambda x: np.log(abs(x)), row)
    log_result.append(log_row)
MS = plt.matshow(log_result)
cbar = plt.colorbar(MS)
titleStr='Wlt analysis with $0\leq j\leq$'+str(J)+' for \n' \
    + r"$\ln\phi(\theta,\theta')=(1-\gamma^2)/(1+\gamma^2-2\gamma\cos(\theta-\theta'))}$"+r'$,\ \gamma=$'+str(g)
plt.title(titleStr)
plt.xlabel('$(j,k)$')
plt.ylabel('$(j\',k\')$')
# plt.show()

titleSave="wlt_analysis-J_"+str(J)+"-gamma_"+str(g)+"-log.pdf"
plt.savefig(titleSave, bbox_inches="tight")
plt.clf()

# ============
# Plot kernel
# ============
# We are using automatic selection of contour levels;
# this is usually not such a good idea, because they don't
# occur on nice boundaries, but we do it here for purposes
# of illustration.
origin = 'lower'
x=np.linspace(-r,r,num=ngrid,endpoint=True)
y=x
X, Y = np.meshgrid(x, y)
PHI=phi(X,Y,g)
# TODO: plot phi as a function of θ-θ'
CS = plt.contourf(X, Y, PHI, 20,
                  alpha=0.8,
                  cmap=plt.cm.bone,
                  origin=origin)
# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(CS)
# cbar.ax.set_ylabel('title of colorbar')

titleStr=r"$\phi(\theta,\theta')=(1-\gamma^2)/(1+\gamma^2-2\gamma\cos(\theta-\theta'))}$"+r'$,\ \gamma=$'+str(g)
plt.title(titleStr)
plt.xlabel(r"$\theta$")
plt.ylabel(r"$\theta'$")
# plt.show()

titleSave="phi"+"-gamma_"+str(g)+".pdf"
plt.savefig(titleSave)

