The file checkLP(t,m,n).m2 is used to check whether QQ[X_1, X_mn]/in(I_t) and QQ[X_1, X_mn]/I_t each have the weak/ strong Lefschetz property, while checkLPinI(t,m,n).m2 is used to check QQ[X_1, X_mn]/in(I_t) only. If we only need to check QQ[X_1, X_mn]/in(I_t) has the weak/ strong Lefschetz property, using checkLPinI(t,m,n).m2 is faster. Moreover, given a homogenous ideal, hasLP.m2 can be used to check if it has the weak/ strong Lefschetz property.

checkLP:

Usage: checkLP(t,m,n)

Inputs: t,m,n, three positive integers

Outputs: whether QQ[X_1, X_mn]/in(I_t) and QQ[X_1, X_mn]/I_t have the weak/ strong Lefschetz property

checkLPinI:

Usage:checkLPinI(t,m,n)

Inputs: t,m,n, three positive integers

Outputs: whether QQ[X_1, X_mn]/in(I_t) has the weak/ strong Lefschetz property

hasLP:

Input: a honogenous ideal inI

Output: whether the quotient ring of inI has the weak/ strong Lefschetz property
