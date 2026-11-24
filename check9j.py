from sympy import sqrt, I, S, simplify, nsimplify, Rational, symbols, re, im, conjugate, expand
from sympy.physics.wigner import wigner_9j, wigner_6j
h={0:-sqrt(3),1:I*sqrt(5)/2}
def Pi(*a):
    p=1
    for x in a: p*=2*x+1
    return sqrt(p)
ReS={1:symbols('R1'),2:symbols('R2')}; ImS={1:symbols('I1'),2:symbols('I2')}
G=symbols('G')  # gamma^{bb0}
for b in (1,2):
  for u in (0,1,2):
    if not (abs(b-1)<=u<=b+1): continue
    tot=0; ok=True
    for g in (0,1,2,3):
      for (l1,l2) in ((1,0),(0,1)):
        lp1=lp2=0; a=0; c=b; v=1; up=b; vp=0
        C=(-1)**(a+l2+lp2-g)*h[l1]*conjugate(h[lp1])*conjugate(h[l2])*h[lp2]*Pi(g,g,b,c,u,v,up,vp)
        C*=wigner_9j(lp1+1,l1+1,g,lp2+1,l2+1,g,b,c,a,prec=None)
        C*=wigner_9j(1,l1,l1+1,1,l2,l2+1,u,v,c,prec=None)
        C*=wigner_9j(1,lp1,lp1+1,1,lp2,lp2+1,up,vp,b,prec=None)
        C=simplify(C)
        if g in (1,2):
          SL1=ReS[g]+I*ImS[g]
          Sval = SL1 if (l1,l2)==(1,0) else conjugate(SL1).subs({conjugate(ReS[g]):ReS[g],conjugate(ImS[g]):ImS[g]})
          tot+=C*Sval
        elif C!=0: ok=False; print('g',g,'nonzero',b,u,C)
    # closed form with two 6j symbols
    dr=0
    for g in (1,2):
      SL1=ReS[g]+I*ImS[g]; SL1c=ReS[g]-I*ImS[g]
      dr+=I*sqrt(5)/2*sqrt((2*b+1)*(2*g+1)*(2*u+1))*wigner_6j(1,1,g,2,1,b)*wigner_6j(2,b,1,u,1,1)*(SL1c+(-1)**(b+u)*SL1)
    print(f'b={b} u={u}  full: {expand(simplify(tot))}   closed-form: {expand(simplify(dr))}')
