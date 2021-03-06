def cathedral(a,b,c,NoPrint=True):
	"""
	Computes the gothic origami corresponding to the cathedral flat surface C(a,b,c) with central cylinder of (even) height c and integral saddle connections a and b as in [McMullen-Mukamel-Wright] \n
	Input:
	
	* a,b,c : lengths (c even)
	
	* NoPrint (opt) : if True, the program does not print nor plots the result
	
	Output:
	
	* s1,s2 : permutations defining the monodromy of the origami
	
	* P : graphic with the tessellation of the cathedral into squares
	\n
	"""
	if mod(c,2) != 0:
		print("The height of the central cylinder must be even!");
		return 0;
	N=0; ## Counter of squares
	l=4*a+2*c; ## length of the central cylinder
	d=c*l+2*a*c+(2*b+c)*c; ## number of squares
	eps=0.5; ## for the display of the number inside the square
	G=SymmetricGroup(d);
	s1=G([()]);
	s2=G([()]);
	
	# Creating s1
	## c central cylinders of width l
	for i in range(l-1):
		for k in range(c):
			s1=G([(k*l+i+1, k*l+i+2)])*s1;
	## b upper and b lower cylinders of width c
	N+=l*c;
	for i in range(c-1):
		for k in range(b):
			s1=G([(c*l+k*c+i+1, c*l+k*c+i+2)])*G([(c*l+b*c+k*c+i+1, c*l+b*c+k*c+i+2)])*s1;
	## c/2 upper and c/2 lower cylinders of width c+2*a
	N+=2*b*c;
	for i in range(2*a+c-1):
		for k in range(c/2):
			s1=G([(N+k*(2*a+c)+i+1, N+k*(2*a+c)+i+2)])*G([(N+(2*a+c)*c/2+k*(2*a+c)+i+1, N+(2*a+c)*c/2+k*(2*a+c)+i+2)])*s1;
	
	# Creating s2
	## central (sub-)cylinders of width c
	for i in range(l):
		for k in range(c-1):
			s2=G([((c-1-k)*l+i+1, (c-2-k)*l+i+1)])*s2;	
	## c left cylinders of width 3*c+2*b
	for i in range(c):
		s2=G([(a+i+1,l*c+i+1)])*s2*G([((c-1)*l+a+i+1,l*c+c*b+i+1)]); ## connector
		for k in range(b-1):
			s2=G([(l*c+k*c+1+i, l*c+(k+1)*c+i+1)])*s2*G([(l*c+b*c+k*c+i+1, l*c+b*c+(k+1)*c+i+1)]);
	for i in range(c):
		s2=G([(l*c+c*(b-1)+i+1,	l*c+2*b*c+(2*a+c)*c/2+(2*a+c)*(c/2-1)+2*a+i+1)])*s2*G([(l*c+2*b*c+(c/2-1)*(c+2*a)+2*a+i+1, l*c+c*(2*b-1)+i+1)]); ## connector
		for k in range(c/2-1):
			s2=G([(l*c+2*b*c+a*c+c*c/2+2*a+(c/2-1-k)*(c+2*a)+i+1,l*c+2*b*c+a*c+c*c/2+2*a+(c/2-1-k-1)*(c+2*a)+i+1)])*s2*G([(l*c+2*b*c+(c/2-1-k)*(c+2*a)+2*a+i+1, l*c+2*b*c+(c/2-1-k-1)*(c+2*a)+2*a+i+1)]);
		## Connecting cylinders
	for i in range(c/2):
		s2=G([(l*c+2*b*c+(2*a+c)*c/2+2*a+i+1,l-c/2+i+1)])*G([(l*c+2*b*c+(2*a+c)*c/2+2*a+c/2+i+1,l-c-2*a+i+1)])*s2;
	## 2*a right cylinders of width 2*c
	for i in range(2*a):
		s2=G([(2*a+c+c/2+i+1, l*c+2*b*c+i+1)])*s2*G([((c-1)*l+2*a+c+c/2+i+1, l*c+2*b*c+(2*a+c)*c/2+i+1)]); ## connector
		for k in range(c/2-1):
			s2=G([(l*c+2*b*c+k*(c+2*a)+i+1,l*c+2*b*c+(k+1)*(c+2*a)+i+1)])*s2*G([(l*c+2*b*c+(2*a+c)*c/2+k*(c+2*a)+i+1, l*c+2*b*c+(2*a+c)*c/2+(k+1)*(c+2*a)+i+1)]);

	if (not(NoPrint)):
	    print("\n", s1, "\n", s2, "\n");
	
	P1=[];
	P2=[];
	P3=[];
	L1=plot([]);
	L2=plot([]);
	L3=plot([]);
	for i in range(l+1):
		for k in range(c+1):
			P1.append([i,k-c/2]);
			if i==0:
				L1+=line([[0,c/2-k],[l,c/2-k]]);
			if k==0:
				L1+=line([[i,-c/2],[i,c/2]]);
	for i in range(c+1):
		for k in range(b+c/2):
			P2.append([a+i,c/2+k+1]);
			P2.append([a+i,-c/2-k-1]);
			if i==0:
				L2+=line([[a,c/2+k+1],[a+c,c/2+k+1]]);
				L2+=line([[a,-c/2-k-1],[a+c,-c/2-k-1]]);
			if k==0:
				L2+=line([[a+i,c/2],[a+i,c+b]]);
				L2+=line([[a+i,-c/2],[a+i,-c-b]]);
	for i in range(2*a+1):
		for k in range(c/2):
			P3.append([2*a+c+c/2+i,c/2+k+1]);
			P3.append([2*a+c+c/2+i,-c/2-k-1]);
			if i==0:
				L3+=line([[2*a+c+c/2,c/2+k+1],[2*a+c+c/2+2*a,c/2+k+1]]);
				L3+=line([[2*a+c+c/2,-c/2-k-1],[2*a+c+c/2+2*a,-c/2-k-1]]);
			if k==0:
				L3+=line([[2*a+c+c/2+i,c/2],[2*a+c+c/2+i,c]]);
				L3+=line([[2*a+c+c/2+i,-c/2],[2*a+c+c/2+i,-c]]);	
	T=[ [i,(-1,0)] for i in [1..d]]; T[0]=[1,(0,c/2-1)];
	Filled=[v[0] for v in T if (v[1][0]!=-1)];
	Unfilled=[v[0] for v in T if (v[1][0]==-1)];
	while Unfilled!=[]: 
		i=Filled[0];
		V=T[i-1][1];	
		if IsInCathedral(a,b,c,V[0]+1,V[1]) & (s1(i) in Unfilled):
			T[s1(i)-1] = [s1(i),(V[0]+1,V[1])];
			Unfilled.remove(s1(i));
			Filled.append(s1(i));
		if IsInCathedral(a,b,c,V[0]-1,V[1]) & ((s1^-1)(i) in Unfilled):
			T[(s1^-1)(i)-1] = [(s1^-1)(i),(V[0]-1,V[1])];
			Unfilled.remove((s1^-1)(i));
			Filled.append((s1^-1)(i));
		if IsInCathedral(a,b,c,V[0],V[1]+1) & (s2(i) in Unfilled):
			T[s2(i)-1] = [s2(i),(V[0],V[1]+1)];
			Unfilled.remove(s2(i));
			Filled.append(s2(i));
		if IsInCathedral(a,b,c,V[0],V[1]-1) & ((s2^-1)(i) in Unfilled):
			T[(s2^-1)(i)-1] = [(s2^-1)(i),(V[0],V[1]-1)];
			Unfilled.remove((s2^-1)(i));
			Filled.append((s2^-1)(i));
		Filled.remove(i);
	txt=sum(plot(text(str(t[0]),(t[1][0]+eps,t[1][1]+eps),fontsize=10, aspect_ratio=1, horizontal_alignment="center", vertical_alignment="center")) for t in T);
	P=point2d(P1+P2+P3,size=20)+L1+L2+L3+txt;
	if (not(NoPrint)):
	    P.show(axes=False, aspect_ratio=1);
	return [s1,s2,P];

def IsInCathedral(a,b,c,X,Y):
	if mod(c,2) != 0:
		print("The height of the central cylinder must be even!");
		return 0;
	l=4*a+2*c; ## length of the central cylinder
	d=c*l+2*a*c+(2*b+c)*c; ## number of squares
	if (X<0) or (X>l-1) or (X<a and (Y>c/2-1 or Y<-c/2)) or ( (X>=a and X<a+c) and (Y>b+c-1 or Y<-b-c)) or ((X>=a+c and X<2*a+c+c/2) and (Y>c/2-1 or Y<-c/2)) or ((X>=2*a+c+c/2 and X<4*a+c+c/2) and (Y>c-1 or Y<-c)) or ((X>=4*a+c+c/2 and X<l) and (Y>c/2-1 or Y<-c/2)):
		return False;
	else:
		return True;


def STable(a,b,c,NoPrint=True):
	"""
    NOT FINISHED
	Computes the S-shaped table origami corresponding to the cathedral flat surface C(a,b,c) with central cylinder of height c and rectangular cylinders with integral saddle connections c and a and width b as in [McMullen] \n
	Input:
	
	* a,b,c : lengths
	
	* NoPrint (opt) : if True, the program does not print nor plots the result
	
	Output:
	
	* s1,s2 : permutations defining the monodromy of the origami
	
	* P : graphic with the tessellation of the S-shaped table into squares (TODO)
	\n
	"""
	N=0; ## Counter of squares
	l=a+c; ## length of the rectangular cylinders
	d=2*l*b+c^2; ## number of squares
	eps=0.5; ## for the display of the number inside the square
	G=SymmetricGroup(d);
	s1=G([()]);
	s2=G([()]);
	
	# Creating s1
	## c central cylinders of width c (central square)
	for i in range(c-1):
		for k in range(c):
			s1=G([(k*c+i+1, k*c+i+2)])*s1;
	## b upper and b lower cylinders of width l
	N+=c*c;
	for i in range(l-1):
		for k in range(b):
			s1=G([(c*c+k*l+i+1, c*c+k*l+i+2)])*G([(c*c+b*l+k*l+i+1, c*c+b*l+k*l+i+2)])*s1;
	
	# Creating s2
	## central (sub-)cylinders of width c
	s2=G([()]);
	for k in range(c):
		for i in range(c):
			s2=G([(k*c+i+1, k*c+c+i+1)])*s2;	
	## c upper cylinders of width b
	for k in range(b-1):
		for i in range(l):
			s2=G([(c*c+k*l+i+1, c*c+k*l+l+i+1)])*s2;	
			s2=G([(c*c+l*b+k*l+i+1, c*c+l*b+k*l+l+i+1)])*s2;	
		## Connecting cylinders
	for i in range(c):
		s2=s2*G([(i+1,c*c+b*l+a+i+1)]);
	return [s1,s2];
