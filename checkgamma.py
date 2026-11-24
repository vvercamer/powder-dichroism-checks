import numpy as np
from sympy.physics.wigner import clebsch_gordan
from functools import lru_cache
@lru_cache(None)
def cg(a,al,b,be,c,ga): return float(clebsch_gordan(a,b,c,al,be,ga))
def sph(v):
    x,y,z=v; return {1:-(x+1j*y)/np.sqrt(2),0:z+0j,-1:(x-1j*y)/np.sqrt(2)}
def couple(A,ra,B,rb,c):
    return {g:sum(cg(ra,a,rb,b2,c,g)*A[a]*B[b2] for a in A for b2 in B if a+b2==g) for g in range(-c,c+1)}
rng=np.random.default_rng(1)
def unit(v): return v/np.linalg.norm(v)
def perp_pol(k):
    a=unit(np.cross(k,rng.normal(size=3))); b=np.cross(k,a)
    z=rng.normal(size=2)+1j*rng.normal(size=2); z/=np.linalg.norm(z)
    return z[0]*a+z[1]*b
for trial in range(3):
    k=unit(rng.normal(size=3)); ks=unit(rng.normal(size=3))
    e=perp_pol(k); es=perp_pol(ks); ec=e.conj(); esc=es.conj()
    Pc=(-1j*np.dot(np.cross(e,ec),k)).real; Pcs=(-1j*np.dot(np.cross(es,esc),ks)).real
    E,Ec,Es,Esc,K=map(sph,(e,ec,es,esc,k))
    kks=np.dot(k,ks); kes2=abs(np.dot(k,es))**2
    def gam(b,u,Out=None):
        In=couple(couple(E,1,Ec,1,u),u,K,1,b)
        if Out is None: Out=couple(Esc,1,Es,1,b)
        return couple(In,b,Out,b,0)[0]
    four=(np.dot(e,esc)*np.dot(np.cross(ec,es),k)+np.dot(e,es)*np.dot(np.cross(ec,esc),k)
          +np.dot(ec,esc)*np.dot(np.cross(e,es),k)+np.dot(ec,es)*np.dot(np.cross(e,esc),k))
    print('trial',trial)
    print(' (1,0)',gam(1,0), Pcs*kks/(3*np.sqrt(2)))
    print(' (1,1)',gam(1,1),0)
    print(' (1,2)',gam(1,2),-Pcs*kks/(3*np.sqrt(10)))
    print(' (2,1)',gam(2,1),Pc/np.sqrt(10)*(1/3-kes2))
    print(' (2,2)',gam(2,2),-1j/(2*np.sqrt(30))*four)
    # average over scattered polarisation: two orthogonal linear pols
    a=unit(np.cross(ks,rng.normal(size=3))); bb=np.cross(ks,a)
    def avg(b,u):
        s=0
        for p in (a,bb):
            P=sph(p); s+=gam(b,u,couple(P,1,P,1,b))
        return s/2
    two=np.dot(e,ks)*np.dot(np.cross(ec,ks),k)+np.dot(ec,ks)*np.dot(np.cross(e,ks),k)
    print(' <1,0>',avg(1,0),' <2,1>',avg(2,1),-Pc/(6*np.sqrt(10))*(1-3*kks**2))
    print(' <2,2>',avg(2,2),1j/(2*np.sqrt(30))*two)
# frame formula
chi,psi,al,be=0.3,0.7,1.1,2.0
k=np.array([0,0,1.]); e=np.array([np.cos(chi)*np.cos(psi)+1j*np.sin(chi)*np.sin(psi),np.cos(chi)*np.sin(psi)-1j*np.sin(chi)*np.cos(psi),0])
ks=np.array([np.sin(al)*np.cos(be),np.sin(al)*np.sin(be),np.cos(al)])
two=np.dot(e,ks)*np.dot(np.cross(e.conj(),ks),k)+np.dot(e.conj(),ks)*np.dot(np.cross(e,ks),k)
print('frame',two, np.cos(2*chi)*np.sin(al)**2*np.sin(2*(be-psi)))
print('e x e* ', np.cross(e,e.conj()), 'i sin2chi=',np.sin(2*chi))
