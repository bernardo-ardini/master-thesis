syms E1 E2 G P dP P0 u1 u2 u3 p1 p2 p3 eta real
syms W(E1,E2,G) H(P,dP) D(P0,P)

W1=diff(W,E1);
W2=diff(W,E2);
W3=diff(W,G);
H1=diff(H,P);
H2=diff(H,dP);
D2=subs(diff(D,P),P0,P);

u=[u1;u2;u3];
p=[p1;p2;p3];

argE1=eta*p(1)/u(3);
argE2=u(2)*u(3);
argG=p(2)/u(3);

W1=formula(subs(W1,[E1,E2,G],[argE1,argE2,argG]));
W2=formula(subs(W2,[E1,E2,G],[argE1,argE2,argG]));
W3=formula(subs(W3,[E1,E2,G],[argE1,argE2,argG]));
H1=formula(subs(H1,[P,dP],[u(3),p(3)]));
H2=formula(subs(H2,[P,dP],[u(3),p(3)]));
D2=formula(D2);

f2=-W2;
f3=(eta*p(1)/u(3)^2)*W1-u(2)*W2+(p(2)/u(3)^2)*W3-H1-D2;
f=[0;f2;f3];
g1=W1*eta/u(3);
g2=W3/u(3);
g3=H2;
g=[g1;g2;g3];

MatA=subs(jacobian(f,u),[p1,p2,p3],[1,0,0]);
MatB=subs(jacobian(f,p),[p1,p2,p3],[1,0,0]);
MatC=subs(jacobian(g,u),[p1,p2,p3],[1,0,0]);
MatD=subs(jacobian(g,p),[p1,p2,p3],[1,0,0]);

W23=formula(subs(diff(W,E2,G),[E1,E2,G],[argE1,argE2,argG]));
W31=formula(subs(diff(W,G,E1),[E1,E2,G],[argE1,argE2,argG]));
W32=formula(subs(diff(W,G,E2),[E1,E2,G],[argE1,argE2,argG]));
H12=formula(subs(diff(H,P,dP),[P,dP],[u3,p3]));
H21=formula(subs(diff(H,dP,P),[P,dP],[u3,p3]));

old=subs([W13,W23,W31,W32,H12,H21],[p1,p2,p3],[1,0,0]);
new=[0,0,0,0,0,0];