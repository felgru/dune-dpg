#!/usr/bin/python3

from __future__ import (absolute_import, print_function)
import argparse
import numpy as np
from scipy.linalg import norm # to compute norm of an array
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

# Scaling function $k_{0,0}$ on the interval [-r,r[
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
    return (x, (xmax-xmin)/(ngrid-1))

# ============
# Main program
# ============
aparser = argparse.ArgumentParser(
        description='plot Henyey–Greenstein kernel')
aparser.add_argument('--no-title', dest='plot_title',
                     action='store_false', default=True,
                     help='plot without title text')
aparser.add_argument('--plot-matrices', dest='plot_matrices',
                     action='store_true', default=False,
                     help='plot matrix representations of the kernel')
aparser.add_argument('--presentation', dest='presentation',
                     action='store_true', default=False,
                     help='generate plot for presentations')
aparser.add_argument('--level', action='store', type=int, default=7,
                     help='maximal wavelet level')
aparser.add_argument('--ngrid', dest='ngrid',
                     action='store', type=int, default=100,
                     help='number of quadrature points')
aparser.add_argument('gamma', action='store', type=float,
                     help='parameter of the Henyey–Greenstein kernel')
args = aparser.parse_args()
if args.presentation:
    args.plot_title = False
ngrid = args.ngrid # Size mesh
J = args.level
r = np.pi
g = args.gamma # forward peak coef

if args.presentation:
    mpl.rc('font',**{'family':'serif','size':20})
mpl.rc('text', usetex=True)

if args.plot_matrices:
    nwlt=0
    result = []
    (y, yquadweight) = quadpoints(0, 0, r)
    # get sf(y)
    sfy=sf(y,r)
    # print(np.sqrt((ymax-ymin)*np.mean(wlty*wlty))) # L2 norm wlt
    # initialize resultx list
    resultx = []
    (x, xquadweight) = quadpoints(0, 0, r)
    # get sf(x)
    sfx=sf(x,r)
    # Values of the kernel
    X, Y = np.meshgrid(x, y)
    PHI=phi(X,Y,g)
    # Product quadweight*kernel*wlt(x)*wlt(y)
    # TODO: in the summed trapezoidal rule, the outer
    #       quadrature points have to be weighted with
    #       factor 1/2.
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
            # TODO: in the summed trapezoidal rule, the outer
            #       quadrature points have to be weighted with
            #       factor 1/2.
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
            # get sf(x)
            sfx=sf(x,r)
            # Values of the kernel
            X, Y = np.meshgrid(x, y)
            PHI=phi(X,Y,g)
            # Product quadweight*kernel*wlt(x)*wlt(y)
            # TODO: in the summed trapezoidal rule, the outer
            #       quadrature points have to be weighted with
            #       factor 1/2.
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
                    # TODO: in the summed trapezoidal rule, the outer
                    #       quadrature points have to be weighted with
                    #       factor 1/2.
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
    MS = plt.matshow(result, cmap=plt.cm.autumn)
    cbar = plt.colorbar(MS)
    if(args.plot_title):
        titleStr = (
            rf'Wlt analysis with $0\leq j\leq {J}$ for \n'
            r'$k(\varphi)=(1-\gamma^2)/'
            r'(1+\gamma^2-2\gamma\cos(\varphi))$'
            rf'$,\ \gamma={g}$'
        )
        plt.title(titleStr)
    plt.xlabel('$(j,k)$')
    plt.ylabel("$(j',k')$")
    # plt.show()

    titleSave = f"wlt_analysis-J_{J}-gamma_{g}.pdf"
    plt.savefig(titleSave, bbox_inches="tight")
    plt.clf()

    # plot logarithmic wavelet representation
    log_result = []
    for row in result:
        log_row = [abs(x) for x in row]
        log_result.append(log_row)
    MS = plt.matshow(log_result, norm=mpl.colors.LogNorm(), cmap=plt.cm.autumn)
    cbar = plt.colorbar(MS)
    if(args.plot_title):
        titleStr = (
            rf'Wlt analysis with $0\leq j\leq {J}$ for \n'
            r'$k(\varphi)=(1-\gamma^2)/'
            r'(1+\gamma^2-2\gamma\cos(\varphi))$'
            rf'$,\ \gamma={g}$'
        )
        plt.title(titleStr)
    plt.xlabel('$(j,k)$')
    plt.ylabel("$(j',k')$")
    # plt.show()

    titleSave=f"wlt_analysis-J_{J}-gamma_{g}-log.pdf"
    plt.savefig(titleSave, bbox_inches="tight")
    plt.clf()

# ============
# Plot kernel
# ============
x=np.linspace(-r,r,num=ngrid,endpoint=True)
# plot phi as a function of θ-θ'
PHI=phi(x, 0, g)
if args.presentation:
    # plot in RWTH blue
    plt.plot(x, PHI, linewidth=2.0, color='#0054AF')
    plt.xlim(-r, r)
    plt.ylim(1e-3, 1e1)
    plt.yscale('log')
else:
    plt.plot(x, PHI)
    plt.xlim(-r, r)

if(args.plot_title):
    titleStr = (
        r'$k(\varphi)=(1-\gamma^2)/'
        r'(1+\gamma^2-2\gamma\cos(\varphi))$'
        rf'$,\ \gamma={g}$'
    )
    plt.title(titleStr)
    plt.ylabel(r"$k(\varphi)$")
plt.xlabel(r"$\varphi$")
# plt.show()

plt.tight_layout()
titleSave = f"phi-gamma_{g}.pdf"
plt.savefig(titleSave)
