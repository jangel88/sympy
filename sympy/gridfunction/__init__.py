from .indexed import IndexedBase, Idx, Indexed
from .array import (MutableDenseNDimArray, ImmutableDenseNDimArray,
    MutableSparseNDimArray, ImmutableSparseNDimArray, NDimArray, tensorproduct,
    tensorcontraction, derive_by_array, permutedims, Array, DenseNDimArray,
    SparseNDimArray,)
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
def Delpx(u,dx):
        return dx*Dpx(u,dx)
def Delmx(u,dx):
        return dx*Dmx(u,dx)
def Apx(u,dx):
	return Rational(1,2)*(u.subs({ix:ix+1})+u.subs({ix:ix}))
def Amx(u,dx):
	return Rational(1,2)*(u.subs({ix:ix})+u.subs({ix:ix-1}))

def Dpy(u,dy):
	return 1/dy*(u.subs({iy:iy+1})-u.subs({iy:iy}))

def Dmy(u,dy):
	return 1/dy*(u.subs({iy:iy})-u.subs({iy:iy-1}))

def D0y(u,dy):
	return 1/(2*dy)*(u.subs({iy:iy+1})-u.subs({iy:iy-1}))
