from .indexed import IndexedBase, Idx, Indexed
#from .array import (MutableDenseNDimArray, ImmutableDenseNDimArray,
#    MutableSparseNDimArray, ImmutableSparseNDimArray, NDimArray, tensorproduct,
#    tensorcontraction, derive_by_array, permutedims, Array, DenseNDimArray,
#    SparseNDimArray,)
from sympy import Rational
from sympy import symbols
ix,iy,iz = symbols('ix iy iz',cld=Idx) 
dx,dy,dz = symbols('dx dy dz')
dr,ds,dt = symbols('dr ds dt')
xi = symbols('xi')
u = IndexedBase('u')
um = IndexedBase('um')
un = IndexedBase('un')
up = IndexedBase('up')

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
	return 1/(2*dr)*(u.subs({ix:ix+1})-u.subs({ix:ix-1}))
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
def Apr(u,dr):
        return Rational(1,2)*(u.subs({ix:ix+1})+u.subs({ix:ix}))
def Amr(u,dr):
        return Rational(1,2)*(u.subs({ix:ix})+u.subs({ix:ix-1}))
def Aps(u,ds):
        return Rational(1,2)*(u.subs({iy:iy+1})+u.subs({iy:iy}))
def Ams(u,ds):
        return Rational(1,2)*(u.subs({iy:iy})+u.subs({iy:iy-1}))

def Dps(u,ds):
	return 1/ds*(u.subs({iy:iy+1})-u.subs({iy:iy}))

def Dms(u,ds):
	return 1/ds*(u.subs({iy:iy})-u.subs({iy:iy-1}))

def D0s(u,ds):
	return 1/(2*ds)*(u.subs({iy:iy+1})-u.subs({iy:iy-1}))

