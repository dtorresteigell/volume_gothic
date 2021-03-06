from surface_dynamics.all import *

load("gothic_cathedral.sage")

N=30; ## ≤120 squares

## Table with the Lyapunov exponents of all Gothic origamis with d≤4*N squares

print("  d    : Class - Approx. of Lyapunov exponents - Index   - Cusps - Sum of Lyapunov exponents\n");
check=0;
for dd in [1..N]:
  Ori=[];
  d=4*dd;
  if check!=0:
    print(" ");
  check=0;
  for cc in divisors(d):
    V=[ [a,b,2*cc] for a in [1..floor((dd/cc-3*cc)/3)] for b in [1..floor(dd/cc-3*cc)] if 2*cc*(6*a+2*b+3*2*cc)==d];
    if V!=[]:
      check=1;
    for v in V:
      [s0,s1,P]=cathedral(v[0],v[1],v[2],True);
      O=Origami(s0,s1);
      Approx=O.lyapunov_exponents_approx();
      G=O.veech_group();
      B=exists(Ori, lambda x : O.is_isomorphic(x));
      if B[0]:
        IsoClass=Ori.index(B[1]);
      else:
        Ori.append(O);
        IsoClass=Ori.index(O);
      print("d=",d," :  ",IsoClass, "  - ", [w.n(digits=5) for w in Approx], " - ",  G.ncusps(), "   - ", G.index(), "  - ", O.sum_of_lyapunov_exponents());