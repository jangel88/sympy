from .indexed import IndexedBase, Idx, Indexed
#from .array import (MutableDenseNDimArray, ImmutableDenseNDimArray,
#    MutableSparseNDimArray, ImmutableSparseNDimArray, NDimArray, tensorproduct,
#    tensorcontraction, derive_by_array, permutedims, Array, DenseNDimArray,
#    SparseNDimArray,)
from sympy import Rational
from sympy import symbols
ix,iy,iz = symbols('ix iy iz',cld=Idx)
dx,dy,dz = symbols('dx dy dz')
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
