-- Hilbert function (Mats Boij)
HF = (I,n) -> (
    A := (ring I)/I;
    m := ideal vars A;
    for i to n list numcols basis ( m^i/m^(i+1) )
    )

--Verify WLP and SLP for R/in(I_t)
checkLP=(t,m,n)->(
    N=m*n;
    R=QQ[x_1..x_N];
    X=transpose genericMatrix(R,x_1,n,m);
    I=minors(t,X);
    inI=ideal leadTerm I;
    d=dim R-codim inI;
    theta=ideal(0);
    for i from 1 to d do (
	Y=random(1,R);
	theta=ideal(Y)+theta;
	);
    J=monomialIdeal leadTerm (inI+theta);
    S=QQ[x_(d+1)..x_N];
    JS=sub (J,S);
    a=max degree numerator reduceHilbert hilbertSeries (JS);
    hf=HF(JS,a);
    print "HF(R/(in(I_t)+theta))=HF(R/(I_t+theta)) is";
    print  hf;
    l=random(1,R);
    c=1;
    cJ=0;
    for s from 1 to a do (
	print ("s =");
	print c;
	use R;
	thetal=ideal(l)^s+theta;
	ILs=monomialIdeal leadTerm (I+thetal);
	JLs=monomialIdeal leadTerm (inI+thetal);
	use S;
	ILsS=sub(ILs,S);
	JLsS=sub(JLs,S);
	H0={1};
	if s>1 then (
	    for j from 1 to s-1 do (
		h=hf#(j);
		H0=H0|{h};
		);
	    );
	for i from 0 to a-s do (
	    h=max{0,hf#(i+s)-hf#i};
	    H0=H0|{h};
	    );
	print "HF(R/(I_t+theta+L^s)) is";
	print  HF(ILsS,a);
	if cJ==0 then (
	    print "HF(R/(in(I_t)+theta+L^s)) = is";
	    print  HF(JLsS,a);
	    if H0!=HF(JLsS,a) then (cJ=c);
		);
	if H0==HF(ILsS,a) then c=c+1 else (
	    if c!=1 then (
		if cJ==1 then (
		    return "R/in(I_t) fails the WLP, R/I_t has the WLP but fails the SLP"
		    ) else (
		    return "Both R/in(I_t) and R/I_t have the WLP but fails the SLP"
		    );
		) else (
		return "Both R/in(I_t) and R/I_t fail the WLP"
		);
	    );
	);
    if a==(c-1) then (
	if cJ==0 then (
	    return "Both R/in(I_t) and R/I_t have the SLP"
	    ) else (
	    if cJ==1 then (
		return "R/in(I_t) fails the WLP, R/I_t has the SLP"
		) else (
		return "R/in(I_t) has the WLP but fails the SLP, R/I_t has the SLP"
		);
	    );	
	);
    )                                                                                                                                                                                                
