clear;
close all;

syms c1 c2 c3 a1 a2 a3 b1 b2 L h reals;

syms v10(x) v20(x);
dv20=diff(v20,x);
dx=dsolve([a1*diff(v10,x,x)+b1*diff(v20,x)==0,a2*diff(v20,x,x)-b1*diff(v10,x)+c1*v20==0],[v10(L)==0,dv20(L)==0]);

syms v1p(x) v2p(x) v3p(x);
dv2p=diff(v2p,x);
dv3p=diff(v3p,x);
sx=dsolve([a1*diff(v1p,x,x)+b1*diff(v2p,x)+b2*diff(v3p,x)==0,a2*diff(v2p,x,x)-b1*diff(v1p,x)+c1*v2p+c2*v3p==0,a3*diff(v3p,x,x)-b2*diff(v1p,x)+c2*v2p+c3*v3p==0],[v1p(0)==0,dv2p(0)==0,dv3p(0)==0]);

syms C1 C2 C3;

matdx=equationsToMatrix(subs([dx.v10;dx.v20;0;diff(dx.v10,x);diff(dx.v20,x)],x,h),[C1,C2]);
matsx=equationsToMatrix(subs([sx.v1p;sx.v2p;sx.v3p;diff(sx.v1p,x);diff(sx.v2p,x)],x,h),[C1,C2,C3]);
Q=[matdx,-matsx];

[~,U]=lu(Q);
de=prod(diag(U));