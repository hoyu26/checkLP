-- Hilbert function (Mats Boij)
HF = (I,n) -> (
    A := (ring I)/I;
    m := ideal vars A;
    for i to n list numcols basis ( m^i/m^(i+1) )
    )

--Verify WLP and SLP for R/in(I_t) and for R/I_t
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
    A=QQ[x_(d+1)..x_N];
    JA=sub (J,A);
    a=max degree numerator reduceHilbert hilbertSeries (JA);
    hf=HF(JA,a);
    l=random(1,R);
    s=1;
    cJ=0;
    for v from 1 to a do (
	use R;
	thetal=ideal(l)^v+theta;
	ILs=monomialIdeal leadTerm (I+thetal);
	JLs=monomialIdeal leadTerm (inI+thetal);
	use A;
	ILsA=sub(ILs,A);
	JLsA=sub(JLs,A);
	H0={1};
	if v>1 then (
	    for j from 1 to v-1 do (
		h=hf#(j);
		H0=H0|{h};
		);
	    );
	for i from 0 to a-v do (
	    h=max{0,hf#(i+s)-hf#i};
	    H0=H0|{h};
	    );
	if cJ==0 then (
	    if H0!=HF(JLsA,a) then (cJ=s);
		);
	if H0==HF(ILsA,a) then s=s+1 else (
	    if s!=1 then (
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
    if a==(s-1) then (
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
