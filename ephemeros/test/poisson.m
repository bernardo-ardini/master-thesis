clear;

addpath('..');

R=3;
C=[0,0];
model=createpde();
circ=[1;C(1);C(2);R];
geometryFromEdges(model,decsg(circ));
generateMesh(model,'Hmax',0.2,'GeometricOrder','linear');
geo=geometry();
geo.readmodel(model);

Omega=geo.whole;
GammaD=Omega.bound.sub(@(x) (abs(x(1,:).^2+x(2,:).^2-R^2)<geo.tol));
GammaN=Omega.bound.sub(@(x) (abs(x(1,:).^2+x(2,:).^2-R^2)<geo.tol)&(x(2,:)<geo.tol));
tic
U=space(Omega,2);
R=space(GammaD,0);
H=space(GammaN,0);
toc
u=field(U);
r=field(R);

alpha=1;

f=field(U);
f.map(@(x) -4*ones(size(x)));
f.eval;

g=field(R);
g.map(@(x) x(1,:).^2+x(2,:).^2);
g.eval;

h=field(H);
h.map(@(x) 2*3*ones(size(x)));
h.eval;

int=alpha*squeeze(sum(reshape(U.shape{2,2},2,U.ref.ndof,1,[]).*reshape(U.shape{2,2},2,1,U.ref.ndof,[]),1));
A=asmbmat(int,2,2,2,2,U,U);
l=asmbvec(f.val{2,1}.*U.shape{2,1},2,2,2,U)+asmbvec(h.val{2,1}.*U.shape{1,1},2,1,1,U);

int=reshape(U.shape{1,1},U.ref.ndof,1,[]).*reshape(R.shape{2,1},1,R.ref.ndof,[]);
B=asmbmat(int,2,1,1,2,U,R);
int=g.val{1,1}.*R.shape{1,1};
n=asmbvec(int,1,1,1,R);

mlt=multi;
mlt.add(U);
mlt.add(R);
mlt.resblk{1}=l;
mlt.resblk{2}=n;
mlt.stiffblk{1,1}=A;
mlt.stiffblk{1,2}=B;
mlt.stiffblk{2,1}=B';

stiff=mlt.stiff;
res=mlt.res;

dof=stiff\res;

u.dof=mlt.comp(dof,1);
r.dof=mlt.comp(dof,2);

u.eval;
figure;
theme(gcf,"light");
u.plot(trisurf=1);