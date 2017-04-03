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
def Dx2(u,dx):
    return D0x(u,dx)
def Dx2mh(u,dx):
    return Dmx(u,dx)
def Dx2ph(u,dx):
    return Dpx(u,dx)
def Dx4(u,dx):
    return Dx2(u,dx)-dx**2*Rational(1,6)*D0x(Dpmx(u,dx),dx)
def Dx4mh(u,dx):
    return Dx2mh(u,dx)-dx**2*Rational(1,24)*Dmx(Dpmx(u,dx),dx)
def Dx4ph(u,dx):
    return Dx2ph(u,dx)-dx**2*Rational(1,24)*Dpx(Dpmx(u,dx),dx)
def Dxx2(u,dx):
        return Dpmx(u,dx)
def Dxx2mh(u,dx):
    return D0x(Dmx(u,dx),dx)
def Dxx2ph(u,dx):
    return D0x(Dpx(u,dx,),dx) 
def Dxx4(u,dx):
        return Dxx2(u,dx)-dx**2*Rational(1,12)*Dpmx2(u,dx) 
def Dxx6(u,dx):
        return Dxx4(u,dx)+dx**4*Rational(1,90)*Dpmx3(u,dx)
def Dxxxx2(u,dx):
        return Dpmx2(u,dx)
def Dxxxx4(u,dx):
        return Dxxxx2(u,dx)-dx**2*Rational(1,6)*Dpmx3(u,dx)
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
def Drr2(u,dr):
        return Dpmr(u,dr)
def Drr4(u,dr):
        return Drr2(u,dr)-dr**2*Rational(1,12)*Dpmr2(u,dr)
def Drr6(u,dr):
        return Drr4(u,dr)+dr**4*Rational(1,90)*Dpmr3(u,dr)
def Drrrr2(u,dr):
        return Dpmr2(u,dr)
def Drrrr4(u,dr):
        return Drrrr2(u,dr)-dr**2*Rational(1,6)*Dpmr3(u,dr)
def Dr2(u,dr):
    return D0r(u,dr)
def Dr4(u,dr):
    return Dr2(u,dr)-Rational(1,6)*dr**2*D0r(Dpmr(u,dr),dr)
def Drrr2(u,dr):
    return D0r(Dpmr(u,dr),dr)

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
def Dss2(u,ds):
        return Dpms(u,ds)
def Dss4(u,ds):
        return Dss2(u,ds)-ds**2*Rational(1,12)*Dpms2(u,ds)
def Dss6(u,ds):
        return Dss4(u,ds)+ds**4*Rational(1,90)*Dpms3(u,ds)
def Dssss2(u,ds):
        return Dpms2(u,ds)
def Dssss4(u,ds):
        return Dssss2(u,ds)-ds**2*Rational(1,6)*Dpms3(u,ds)
def Ds2(u,ds):
    return D0s(u,ds)
def Ds4(u,ds):
    return Ds2(u,ds)-Rational(1,6)*ds**2*D0s(Dpms(u,ds),ds)




def Dr2mh(u,dr):
    return Dmr(u,dr)
def Dr2ph(u,dr):
    return Dpr(u,dr)
def Dr4mh(u,dr):
    return Dr2mh(u,dr)-Rational(1,24)*dr**2*Dmr(Dpmr(u,dr),dr)
def Dr4ph(u,dr):
    return Dr2ph(u,dr)-Rational(1,24)*dr**2*Dpr(Dpmr(u,dr),dr) 
def Drr2mh(u,dr):
    return D0r(Dmr(u,dr),dr)
def Drr4mh(u,dr):
    return Drr2mh(u,dr)
def Drrr2mh(u,dr):
    return Dmr(Dpmr(u,dr),dr)
def Drrr2ph(u,dr):
    return Dpr(Dpmr(u,dr),dr)
def Drrr4mh(u,dr): 
    return Drrr2mh(u,dr)-dr**2*Rational(1,8)*Dmr(Dpmr2(u,dr),dr)
def Drrr4ph(u,dr):
    return Drrr2ph(u,dr)-dr**2*Rational(1,8)*Dpr(Dpmr2(u,dr),dr)


def Ds2mh(u,ds):
    return Dms(u,ds)
def Ds2ph(u,ds):
    return Dps(u,ds)
def Ds4mh(u,ds):
    return Ds2mh(u,ds)-Rational(1,24)*ds**2*Dms(Dpms(u,ds),ds)
def Ds4ph(u,ds):
    return Ds2ph(u,ds)-Rational(1,24)*ds**2*Dps(Dpms(u,ds),ds) 
def Dss2mh(u,ds):
    return D0s(Dms(u,ds),ds)
def Dsss2mh(u,ds):
    return Dms(Dpms(u,ds),ds)
def Dsss2ph(u,ds):
    return Dps(Dpms(u,ds),ds)
def Dsss4mh(u,ds): 
    return Dsss2mh(u,ds)-ds**2*Rational(1,8)*Dms(Dpms2(u,ds),ds)
def Dsss4ph(u,ds):
    return Dsss2ph(u,ds)-ds**2*Rational(1,8)*Dps(Dpms2(u,ds),ds)





def Delmr(u,dr):
    return (u.subs({ix:ix})-u.subs({ix:ix-1}))
def Delpr(u,dr):
    return (u.subs({ix:ix+1})-u.subs({ix:ix}))
def Delms(u,ds):
    return (u.subs({iy:iy})-u.subs({iy:iy-1}))
def Delps(u,dr):
    return (u.subs({iy:iy+1})-u.subs({iy:iy}))




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
def Mrs2d2h(udot):
    return (1/J[ix,iy]*Dpr(Rational(1,2)*Amr2h(J[ix,iy]*c*sqrt(g11[ix,iy]),dr)*(Delmr(udot,dr)-dr*D0r(Amr2h(udot,dr),dr)),dr)
           +1/J[ix,iy]*Dps(Rational(1,2)*Ams2h(J[ix,iy]*c*sqrt(g22[ix,iy]),ds)*(Delms(udot,ds)-ds*D0s(Ams2h(udot,ds),ds)),ds)) 
def Mrs2d4h(udot):
    return (1/J[ix,iy]*Dpr(Rational(1,2)*Amr4h(J[ix,iy]*c*sqrt(g11[ix,iy]),dr)*(Delmr(udot,dr)-dr*Dr4(Amr2h(udot,dr),dr))
            +Rational(1,12)*dr**2*Drr4(Delmr(u,dr),dr),dr)
           +1/J[ix,iy]*Dps(Rational(1,2)*Ams4h(J[ix,iy]*c*sqrt(g22[ix,iy]),ds)*(Delms(udot,ds)-ds*Ds4(Ams2h(udot,ds),ds))
            +Rational(1,12)*ds**2*Dss4(Delms(u,ds),ds),ds))
def Lrs2d2h(u):
    return (1/J[ix,iy]*(Dpr(Amr2h(J[ix,iy]*c**2*g11[ix,iy],dr)*Dmr(u,dr)
                                +Amr2h(J[ix,iy]*c**2*g12[ix,iy],dr)*Amr2h(D0s(u,ds),dr),dr)
                       +Dps(Ams2h(J[ix,iy]*c**2*g21[ix,iy],ds)*Ams2h(D0r(u,dr),ds)
                                +Ams2h(J[ix,iy]*c**2*g22[ix,iy],ds)*Dms(u,ds),ds)))
def Lrs2d4h(u):
    return (1/J[ix,iy]*(Dpr(Amr4h(J[ix,iy]*c**2*g11[ix,iy],dr)*Dr4mh(u,dr)
                           +Amr4h(J[ix,iy]*c**2*g12[ix,iy],dr)*Amr4h(Ds4(u,ds),dr)
                                -dr**2*Rational(1,24)*(Drr2mh(J[ix,iy]*c**2*g11[ix,iy],dr)*Dr4mh(u,dr)
                                                       +Amr4h(J[ix,iy]*c**2*g11[ix,iy],dr)*Drrr2mh(u,dr)
                                                       +Drr2mh(J[ix,iy]*c**2*g12[ix,iy],dr)*Amr4h(Ds4(u,ds),dr)
                                                       +Amr4h(J[ix,iy]*c**2*g12[ix,iy],dr)*Drr2mh(Ds4(u,ds),dr)),dr)
                       +Dps(Ams4h(J[ix,iy]*c**2*g21[ix,iy],ds)*Ams4h(Dr4(u,dr),ds)
                           +Ams4h(J[ix,iy]*c**2*g22[ix,iy],ds)*Ds4mh(u,ds)
                                -ds**2*Rational(1,24)*(Dss2mh(J[ix,iy]*c**2*g21[ix,iy],ds)*Ams4h(Dr4(u,dr),ds)
                                                       +Ams4h(J[ix,iy]*c**2*g21[ix,iy],ds)*Dss2mh(Dr4(u,dr),ds)
                                                       +Dss2mh(J[ix,iy]*c**2*g22[ix,iy],ds)*Ds4mh(u,ds)
                                                       +Ams4h(J[ix,iy]*c**2*g22[ix,iy],ds)*Dsss2mh(u,ds)),ds) 
                                ))
