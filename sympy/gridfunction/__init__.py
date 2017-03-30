from .indexed import IndexedBase, Idx, Indexed
#from .array import (MutableDenseNDimArray, ImmutableDenseNDimArray,
#    MutableSparseNDimArray, ImmutableSparseNDimArray, NDimArray, tensorproduct,
#    tensorcontraction, derive_by_array, permutedims, Array, DenseNDimArray,
#    SparseNDimArray,)
from sympy import Rational
from sympy import symbols,sqrt
ix,iy,iz = symbols('ix iy iz',cld=Idx) 
dx,dy,dz = symbols('dx dy dz')
dr,ds,dt = symbols('dr ds dt')
xi = symbols('xi')
J = IndexedBase('J')
g11 = IndexedBase('g11')
g12 = IndexedBase('g12')
g21 = IndexedBase('g21')
g22 = IndexedBase('g22')
c = symbols('c')
u = IndexedBase('u')
um = IndexedBase('um')
un = IndexedBase('un')
up = IndexedBase('up')
udot = IndexedBase('udot')

def Dpx(u,dx):
	return 1/dx*(u.subs({ix:ix+1})-u.subs({ix:ix}))

def Dmx(u,dx):
	return 1/dx*(u.subs({ix:ix})-u.subs({ix:ix-1}))

def D0x(u,dx):
	return 1/(2*dx)*(u.subs({ix:ix+1})-u.subs({ix:ix-1}))
def Dpmx(u,dx):
        return Dpx(Dmx(u,dx),dx)
def Dpmx2(u,dx):
        return Dpmx(Dpmx(u,dx),dx) 
def Dpmx3(u,dx):
        return Dpmx(Dpmx2(u,dx),dx) 
def Dpmx4(u,dx):
        return Dpmx(Dpmx3(u,dx),dx) 
def Dxx2h(u,dx):
        return Dpmx(u,dx)
def Dxx4h(u,dx):
        return Dxx2h(u,dx)-dx**2*Rational(1,12)*Dpmx2(u,dx)
def Dxx6h(u,dx):
        return Dxx4h(u,dx)+dx**4*Rational(1,90)*Dpmx3(u,dx)
def Dxxxx2h(u,dx):
        return Dpmx2(u,dx)
def Dxxxx4h(u,dx):
        return Dxxxx2h(u,dx)-dx**2*Rational(1,6)*Dpmx3(u,dx)
def Dpy(u,dy):
	return 1/dy*(u.subs({iy:iy+1})-u.subs({iy:iy}))

def Dmy(u,dy):
	return 1/dy*(u.subs({iy:iy})-u.subs({iy:iy-1}))

def D0y(u,dy):
	return 1/(2*dy)*(u.subs({iy:iy+1})-u.subs({iy:iy-1}))

def Dpr(u,dr):
	return 1/dr*(u.subs({ix:ix+1})-u.subs({ix:ix}))

def Dmr(u,dr):
	return 1/dr*(u.subs({ix:ix})-u.subs({ix:ix-1}))

def D0r(u,dr):
	return Rational(1,2)*1/dr*(u.subs({ix:ix+1})-u.subs({ix:ix-1}))
def Dpmr(u,dr):
        return Dpr(Dmr(u,dr),dr)
def Dpmr2(u,dr):
        return Dpmr(Dpmr(u,dr),dr) 
def Dpmr3(u,dr):
        return Dpmr(Dpmr2(u,dr),dr) 
def Dpmr4(u,dr):
        return Dpmr(Dpmr3(u,dr),dr) 
def Drr2h(u,dr):
        return Dpmr(u,dr)
def Drr4h(u,dr):
        return Drr2h(u,dr)-dr**2*Rational(1,12)*Dpmr2(u,dr)
def Drr6h(u,dr):
        return Drr4h(u,dr)+dr**4*Rational(1,90)*Dpmr3(u,dr)
def Drrrr2h(u,dr):
        return Dpmr2(u,dr)
def Drrrr4h(u,dr):
        return Drrrr2h(u,dr)-dr**2*Rational(1,6)*Dpmr3(u,dr)
def Apr2h(u,dr):
        return Rational(1,2)*(u.subs({ix:ix+1})+u.subs({ix:ix}))
def Amr2h(u,dr):
        return Rational(1,2)*(u.subs({ix:ix})+u.subs({ix:ix-1}))
def Aps2h(u,ds):
        return Rational(1,2)*(u.subs({iy:iy+1})+u.subs({iy:iy}))
def Ams2h(u,ds):
        return Rational(1,2)*(u.subs({iy:iy})+u.subs({iy:iy-1}))
def Apr4h(u,dr):
        return Rational(1,16)*(-u.subs({ix:ix+2})+9*u.subs({ix:ix+1})+9*u.subs({ix:ix})-u.subs({ix:ix-1}))
def Amr4h(u,dr):
        return Rational(1,16)*(-u.subs({ix:ix+1})+9*u.subs({ix:ix})+9*u.subs({ix:ix-1})-u.subs({ix:ix-2}))
def Aps4h(u,ds):
        return Rational(1,16)*(-u.subs({iy:iy+2})+9*u.subs({iy:iy+1})+9*u.subs({iy:iy})-u.subs({iy:iy-1}))
def Ams4h(u,ds):
        return Rational(1,16)*(-u.subs({iy:iy+1})+9*u.subs({iy:iy})+9*u.subs({iy:iy-1})-u.subs({iy:iy-2}))

def Dps(u,ds):
	return 1/ds*(u.subs({iy:iy+1})-u.subs({iy:iy}))

def Dms(u,ds):
	return 1/ds*(u.subs({iy:iy})-u.subs({iy:iy-1}))

def D0s(u,ds):
	return Rational(1,2)*1/ds*(u.subs({iy:iy+1})-u.subs({iy:iy-1}))

def Dpms(u,ds):
        return Dps(Dms(u,ds),ds)
def Dpms2(u,ds):
        return Dpms(Dpms(u,ds),ds) 
def Dpms3(u,ds):
        return Dpms(Dpms2(u,ds),ds) 
def Dpms4(u,ds):
        return Dpms(Dpms3(u,ds),ds) 
def Dss2h(u,ds):
        return Dpms(u,ds)
def Dss4h(u,ds):
        return Dss2h(u,ds)-ds**2*Rational(1,12)*Dpms2(u,ds)
def Dss6h(u,ds):
        return Dss4h(u,ds)+ds**4*Rational(1,90)*Dpms3(u,ds)
def Dssss2h(u,ds):
        return Dpms2(u,ds)
def Dssss4h(u,ds):
        return Dssss2h(u,ds)-ds**2*Rational(1,6)*Dpms3(u,ds)
def Delmr(u,dr):
    return (u.subs({ix:ix})-u.subs({ix:ix-1}))
def Delpr(u,dr):
    return (u.subs({ix:ix+1})-u.subs({ix:ix}))
def Delms(u,ds):
    return (u.subs({iy:iy})-u.subs({iy:iy-1}))
def Delps(u,dr):
    return (u.subs({iy:iy+1})-u.subs({iy:iy}))
def Mrs2d2h(udot):
    return (1/J[ix,iy]*Dpr(Rational(1,2)*Amr2h(J[ix,iy]*c*sqrt(g11[ix,iy]),dr)*(Delmr(udot,dr)-dr*D0r(Amr2h(udot,dr),dr)),dr)
           +1/J[ix,iy]*Dps(Rational(1,2)*Ams2h(J[ix,iy]*c*sqrt(g22[ix,iy]),ds)*(Delms(udot,ds)-ds*D0s(Ams2h(udot,ds),ds)),ds))
def Lrs2d2h(u):
    return (1/J[ix,iy]*(Dpr(Amr2h(J[ix,iy]*c**2*g11[ix,iy],dr)*Dmr(u,dr)
                                +Amr2h(J[ix,iy]*c**2*g12[ix,iy],dr)*Amr2h(D0s(u,ds),dr),dr)
                       +Dps(Ams2h(J[ix,iy]*c**2*g21[ix,iy],ds)*Ams2h(D0r(u,dr),ds)
                                +Ams2h(J[ix,iy]*c**2*g22[ix,iy],ds)*Dms(u,ds),ds)))
def Lrs2d4h(u):
    return (1/J[ix,iy]*(Dpr(Amr4h(J[ix,iy]*c**2*g11[ix,iy],dr)*Dmr(u-dr**2*Rational(1,24)*Dpmr(u,dr),dr)
                                +Amr4h(J[ix,iy]*c**2*g12[ix,iy],dr)*Amr4h(D0s(u-ds**2*Rational(1,6)*Dpms(u,ds),ds),dr)
                                -dr**2*Rational(1,24)*(D0r(Dmr(J[ix,iy]*c**2*g11[ix,iy],dr),dr)*Dmr(u-dr**2*Rational(1,24)*Dpmr(u,dr),dr)
                                                       +Amr4h(J[ix,iy]*c**2*g11[ix,iy],dr)*Dpmr(Dmr(u,dr),dr)
                                                       +D0r(Dmr(J[ix,iy]*c**2*g12[ix,iy],dr),dr)*Amr4h(D0s(u-ds**2*Rational(1,6)*Dpms(u,ds),ds),dr)
                                                       +Amr4h(J[ix,iy]*c**2*g12[ix,iy],dr)*D0r(Dmr(D0s(u-ds**2*Rational(1,6)*Dpms(u,ds),ds),dr),dr)),dr)
                               +Dps(Ams4h(J[ix,iy]*c**2*g21[ix,iy],ds)*Dms(u-ds**2*Rational(1,24)*Dpms(u,ds),ds)
                                +Ams4h(J[ix,iy]*c**2*g22[ix,iy],ds)*Ams4h(D0r(u-dr**2*Rational(1,6)*Dpmr(u,dr),dr),ds)
                                -ds**2*Rational(1,24)*(D0s(Dms(J[ix,iy]*c**2*g21[ix,iy],ds),ds)*Dms(u-ds**2*Rational(1,24)*Dpms(u,ds),ds)
                                                       +Ams4h(J[ix,iy]*c**2*g21[ix,iy],ds)*Dpms(Dms(u,ds),ds)
                                                       +D0s(Dms(J[ix,iy]*c**2*g22[ix,iy],ds),ds)*Ams4h(D0r(u-dr**2*Rational(1,6)*Dpmr(u,dr),dr),ds)
                                                       +Ams4h(J[ix,iy]*c**2*g22[ix,iy],ds)*D0s(Dms(D0r(u-dr**2*Rational(1,6)*Dpmr(u,dr),dr),ds),ds)),ds)
                                ))
