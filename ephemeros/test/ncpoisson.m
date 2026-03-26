clear;

addpath('..');

% topol=[1,2,3;1,3,4];
% coord=[0,0;1,0;1,1;0,1];
% geo=geometry();
% geo.init(topol,coord,gord=0);
% figure;
% theme(gcf,"light");
% geo.plot(label=[1,2,3]);

L=1;
H=1;
model=createpde();
rect=[3;4;0;1;1;0;0;0;H;H];
geometryFromEdges(model,decsg(rect));
generateMesh(model,'Hmax',0.05,'GeometricOrder','linear');
geo=geometry();
tic
geo.readmodel(model);
toc
%figure
%theme(gcf,"light");
%geo.plot(label=[1,2,3]);

Omega=geo.whole;
GammaD=Omega.bound.sub(@(x) (abs(x(2,:))<geo.tol)+(abs(x(1,:))<geo.tol)+(abs(x(1,:)-1)<geo.tol));
GammaN=Omega.bound.sub(@(x) (abs(x(2,:)-1)<geo.tol));
Gamma=Omega.skel;
U=space(Omega,-1);
R=space(GammaD,0);
M=space(Gamma,0);
H=space(GammaN,0);
u=field(U);
r=field(R);
m=field(M);

alpha=1;

f=field(U);
f.map(@(x) -4*ones(size(x)));
f.eval;

g=field(R);
g.map(@(x) (x(1,:)-1/2).^2+(x(2,:)-1/2).^2);
g.eval;

h=field(H);
h.map(@(x) 2*(x(2,:)-1/2));
h.eval;

int=alpha*squeeze(sum(reshape(U.shape{2,2},2,U.ref.ndof,1,[]).*reshape(U.shape{2,2},2,1,U.ref.ndof,[]),1));
A=asmbmat(int,2,2,2,2,U,U);
l=asmbvec(f.val{2,1}.*U.shape{2,1},2,2,2,U)+asmbvec(h.val{2,1}.*U.shape{1,1},2,1,1,U);

int=reshape(U.shape{1,1},U.ref.ndof,1,[]).*reshape(R.shape{2,1},1,R.ref.ndof,[]);
B=asmbmat(int,2,1,1,2,U,R);
int=g.val{1,1}.*R.shape{1,1};
n=asmbvec(int,1,1,1,R);

int=reshape(U.shape{1,1},U.ref.ndof,1,[]).*reshape(M.shape{2,1},1,M.ref.ndof,[]);
int=int.*reshape(geo.gauss{2}.trace{1}.sign,1,1,[]);
C=asmbmat(int,2,1,1,2,U,M);

mlt=multi;
mlt.add(U);
mlt.add(R);
mlt.add(M);
mlt.resblk{1}=l;
mlt.resblk{2}=n;
mlt.stiffblk{1,1}=A;
mlt.stiffblk{1,2}=B;
mlt.stiffblk{2,1}=B';
mlt.stiffblk{1,3}=C;
mlt.stiffblk{3,1}=C';

stiff=mlt.stiff;
res=mlt.res;

dof=stiff\res;

u.dof=mlt.comp(dof,1);
r.dof=mlt.comp(dof,2);
m.dof=mlt.comp(dof,3);

u.eval;
figure;
theme(gcf,"light");
u.plot(trisurf=1);