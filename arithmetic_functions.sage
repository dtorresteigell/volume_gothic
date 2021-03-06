"""
Arithmetic functions from 
  D. Zagier "On the values at negative integers of the zeta-function of a real quadratic ﬁeld"
  M. Möller, D. Torres-Teigell "Euler characteristics of Gothic Teichmüller curves"

Application to the computation of the Masur-Veech volume of the gothic locus and the Prym loci.
"""

def Q(r,n):
  if mod(r,2)==0 and (2^(r-2)).divides(n) and mod(n//(2^(r-2)),4)==1:
    return 2^(r//2)*(-1)^(((n//(2^(r-2)))-1)/4);
  elif mod(r,2)==1 and (2^(r-1)).divides(n):
    m=n//(2^(r-1));
    return 2^((r-1)/2)*(-1)^(m*(m-1)/2);
  else:
    return 0;


def gamma(c,n):
  l=c.squarefree_part();
  d=round(sqrt(c/l));
  if mod(c,2)==0:
    return(Q(c.valuation(2),n)*gamma(c//(2^(c.valuation(2))),n));
  elif d.divides(n)==0:
    return(0);
  else:
    return sum([t*moebius(d/t)*kronecker(n//(t^2),l) for t in divisors(gcd(d,n//d))]);


def gamma3(r,d):
  if r==0 or (r==1 and d%3!=0):
    return 1;
  elif r%2==0 and (d^2)%(3^r)==0:
    return 3^(r//2-1)*2;
  elif r%2==1 and valuation(d^2,3)==r-1:
    return 3^((r-1)//2);
  else:
    return 0;


def gamma2(r,d):
  if r==0:
    return 1;
  elif r%2==0 and valuation(d^2,2)==r-2:
    return 2^(r//2);
  elif r%2==1 and valuation(d^2,2)>=r-1:
    return 2^((r-1)//2);
  else:
    return 0;


def gammap(p,r,d):
  if r==0 or (r==1 and d%p!=0):
    return 1;
  elif r%2==0 and (d^2)%(p^r)==0:
    return p^(r//2-1)*(p-1);
  elif r%2==1 and valuation(d^2,p)==r-1:
    return p^((r-1)//2);
  else:
    return 0;

"""
TEST:
for p in primes(100):
  if p==2: 
    for r in [0..50]:
      for d in [1..100]:
        if gamma(p^r,d^2) -  gamma2(r,d) != 0:
          print p, " , ",d, " , ", r, " :  ", gamma(p^r,d^2), " - " , gamma2(r,d);
  else:
    for r in [0..50]:
      for d in [1..100]:
        if gamma(p^r,d^2) -  gammap(p,r,d) != 0:
          print p, " , ",d, " , ", r, " :  ", gamma(p^r,d^2), " - " , gammap(p,r,d);
"""



def EulerProd(p,n,k=1,S='s',prec=10):
  if S=='s':
    s=var('s');
    S=eval(S);
  return 1+sum([gcd(2*k,p^j)^2*gamma(p^j,n)/p^(j*S) for j in [1..prec]]);


def P(p,D,k=1,S='s'):
  if S=='s':
    s=var('s');
    S=eval(S);
  if p==2:
    if p.divides(D) and D.valuation(p)%2==0:
      return 1 + sum([gcd(2*k,p^(2*j+1))^2*p^(j)/p^((2*j+1)*S) for j in [0..(D.valuation(p))//2] ]) + gcd(2*k,p^(2*((D.valuation(p)//2) +1)))^2*p^((D.valuation(p)//2) +1)/p^(((D.valuation(p))+2)*S)
    elif p.divides(D)==0:
      return 1 + gcd(2*k,p)^2*1/p^S  + gcd(2*k,p^2)^2*2/p^(2*S);
  else:
    if p.divides(D) and D.valuation(p)%2==0:
      a=(p-1)/p^(2*S); r=p^(1-2*S);
      return 1 + gcd(2*k,p)^2*a*(1-r^(D.valuation(p)//2))/(1-r) + gcd(2*k,p)^2*p^(D.valuation(p)//2)/p^(((D.valuation(p))+1)*S);
      # return 1 + sum([p^(k-1)*(p-1)/p^(2*k*S) for k in [1..(D.valuation(p))//2] ]) + p^(D.valuation(p)//2)/p^(((D.valuation(p))+1)*S)
    elif p.divides(D)==0:
      return 1 + gcd(2*k,p)^2*1/p^S;


def Palt(p,D,k=1,S='s'):
  if S=='s':
    s=var('s');
    S=eval(S);
  if p==2:
    if p.divides(D) and k==1 and D.valuation(p)%2==0:
      a=1/p^(S-2); r=p^(1-2*S); t=D.valuation(p)//2;
      return 1 + gcd(2*k,p)^2*(p^(3*t+1)+p^2 + 1)/(p^(3*t+2)*(p^2 + p +1));
      #return 1 + (p^(3*t+2)-p^(3*t+1)-p^2 +p +p^3 -1)/(p^(3*t+2)*(p^3-1));
    if p.divides(D) and D.valuation(p)%2==0:
      a=1/p^(S-2); r=p^(1-2*S); t=D.valuation(p)//2;
      return 1 + gcd(2*k,p)^2*(p^(3*t+1)+p^2 + 1)/(p^(3*t+2)*(p^2 + p +1));
      #return 1 + (p^(3*t+2)-p^(3*t+1)-p^2 +p +p^3 -1)/(p^(3*t+2)*(p^3-1));
    elif p.divides(D)==0:
      return 1 + gcd(2*k,p)^2*1/p^S;
  else:
    if p.divides(D) and D.valuation(p)%2==0:
      a=(p-1)/p^(2*S); r=p^(1-2*S); t=D.valuation(p)//2;
      return 1 + gcd(2*k,p)^2*(p^(3*t+1)+p^2 + 1)/(p^(3*t+2)*(p^2 + p +1));
      #return 1 + (p^(3*t+2)-p^(3*t+1)-p^2 +p +p^3 -1)/(p^(3*t+2)*(p^3-1));
    elif p.divides(D)==0:
      return 1 + gcd(2*k,p)^2*1/p^S;

def e(n,k=1,M=10000,prec=10,Coef=1,Round=0):
  """
  If Coef, it returns the coefficient \bar e(n,k) as in "Euler characteristics of Gothic Teichmüller curves"
  Otherwise, it avoids the D^3/2 factor
  """
  Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2)/k^2;
  aux=Const*prod([P(p,n,k,2) for p in primes(2,M) ]); # == Const*n^(3/2)*prod([EulerProd(p,n,k,2,prec) for p in primes(2,M) ]);
  if Coef:
    aux=aux*n^(3/2);
  if Round:
    return round(aux.n());
  else:
    return aux.n();


def e6(D,M=10000, Round=1):
  Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2)/6^2;
  P1=prod([P(p,D,k=1,S=2) for p in primes(2,M)]);
  P2=prod([P(p,D//2^D.valuation(2),k=1,S=2) for p in primes(2,M)]);
  P3=prod([P(p,D//3^D.valuation(3),k=1,S=2) for p in primes(2,M)]);
  P6=prod([P(p,D//(2^D.valuation(2)*3^D.valuation(3)),k=1,S=2) for p in primes(2,M)]);
  #print round(P1), "\n" , round(P2), "\n" , round(P3), "\n" , round(P6), "\n";
  if Round:
    return round(Const*D^(3/2)*(9*4*P1-9*12/5*P2-36/5*4*P3+36/5*12/5*P6).n());
  else:
    return Const*D^(3/2)*(9*4*P1-9*12/5*P2-36/5*4*P3+36/5*12/5*P6);


def ee6(D,M=10000, RR=1):
  #Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2); # Constant of e(n,k) for k=1;
  P1=e(D,1,M,Coef=0,Round=0);
  P2=e(D//2^D.valuation(2),1,M,Coef=0,Round=0);
  P3=e(D//3^D.valuation(3),1,M,Coef=0,Round=0);
  P6=e(D//(2^D.valuation(2)*3^D.valuation(3)),1,M,Coef=0,Round=0);
  #print P1, "\n" , P2, "\n" , P3, "\n" , P6, "\n";
  if RR:
    return round(D^(3/2)*(P1 - 3/5*P2 - 4/5*P3 + 12/25*P6).n());
  else:
    return D^(3/2)*(P1-3/5*P2-4/5*P3+12/25*P6);



# e(d²,6) = d³ · ( e(d²,1) - 3/5 * e(d_2²,1) - 4/5 * e(d_3²,1) + 12/25 * e(d_6²,1) ) , where e(D,1) lacks the D^(3/2) factor
"""
for D in [1000000..1050000]:
  if is_square(D) and (D%24) in [0,1,4,9,12,16]: 
    print D, "  :  ", e6(D) - e(D,6,Coef=1,Round=1);
"""


def esum2(D,M=10000):
  Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2)/6^2;
  return (Const*D^(3/2)*prod([P(p,D//2^D.valuation(2),k=1,S=2) for p in primes(2,M)])).n();


def esum3(D,M=10000):
  Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2)/6^2;
  return (Const*D^(3/2)*prod([P(p,D//3^D.valuation(3),k=1,S=2) for p in primes(2,M)])).n();


def esum6(D,M=10000):
  Const=pi^(1/2)*zeta(2)/2^4/gamma__exact(5/2)/6^2;
  return (Const*D^(3/2)*prod([P(p,D//(2^D.valuation(2)*3^D.valuation(3)),k=1,S=2) for p in primes(2,M)])).n();




def Comp(D):
  if D%12==0:
    return 1;
  elif D%24==1:
    return 2;
  elif D%24==9:
    return 4/3;
  else:
    return 3/2;


def em(d):
  return sum([moebius(d//f)*round(f^3*e(f^2,1)) for f in divisors(d) ])

def e2m(d):
  return sum([moebius(d//f)*round(f^3*e((f//2^f.valuation(2))^2,1)) for f in divisors(d) ])

def e3m(d):
  return sum([moebius(d//f)*round(f^3*e((f//3^f.valuation(3))^2,1)) for f in divisors(d) ])

def e6m(d):
  return sum([moebius(d//f)*round(f^3*e((f//(3^f.valuation(3)*2^f.valuation(2)))^2,1)) for f in divisors(d) ])


def einv(d,k=1):
  return sum([moebius(f)*e((d/f)^2,k) for f in divisors(d)]);

def minv(f,d):
  return sum([moebius(h)*f(d/h) for h in divisors(d)]);



""" 
Proof of Lemma 5.4

for k in [6..20]:
  if k.is_squarefree():
    for d in [5..15]:
      print "d^2=",d^2, ", k=", k, "  -  ", abs(e(d^2,1,5000,Coef=0)+ sum([moebius(m)*prod([(p^2-1)/(p^2+1) for p in prime_divisors(m)])*e(div(d,m)^2,1,5000,Coef=0) for m in divisors(k) if m!=1]) - e(d^2,k,50000, Coef=0));
"""


def test(d,k):
  S1=sum([moebius(d//m)*m^3*e(div(m,k)^2,1,Coef=0) for m in divisors(d) ]);
  S2=prod([p^(3*d.valuation(p)-3)*(p^3-1) for p in prime_divisors(k) if p.divides(d)])*sum([moebius(div(d,k)//m)*m^3*e(m^2,1,Coef=0) for m in divisors(div(d,k)) ]);
  #print S1, " - ", S2;
  print  "\n d=", d, "; k=", k, "  -  ", round(S1/S2);
