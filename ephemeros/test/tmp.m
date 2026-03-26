A=10*rand(2,2);

disp(norm(A));

[dG,d2G]=pl.mat.diss(A);

epsilon=1e-6;
H=epsilon*rand(2,2);
dG1=pl.mat.diss(A+H);

format shortE

ddG=dG1-dG;
ddGex=tensorprod(d2G,H,[3,4],[1,2]);

disp([ddG,ddGex,ddGex-ddG]);
