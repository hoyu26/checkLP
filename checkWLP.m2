checkWLP=(I)->(
theta=ideal();
d=dim ((ring I)/I);
for i from 1 to d do theta=ideal(random(1,ring I))+theta;
J=ideal leadTerm (I+theta);
A=QQ[drop (flatten entries vars ring I, d)];
JA=trim sub (J,A);
a=max degree numerator reduceHilbert hilbertSeries (JA);
hf=hilF(JA,a);
print ("HF=", hf);
l=random(1,ring I);
use ring I;
thetal=ideal(l)+theta;
ILs=ideal leadTerm (I+thetal);
ILsA=sub(ILs,A);
H1=hilF(ILsA,a);
H0={1};
for i from 0 to a-1 do (
	h=max{0,hf#(i+1)-hf#i};
	H0=H0|{h};
	); 
print ("s=",1);
print ("max{0, nabla^s HF}=", H0);
print ("the Hilbert function of the quotient by L^s is", H1);
if H0!=H1 then print "The chosen linear form is not a LE for this ideal." else print "This ideal has the WLP." ;
)
