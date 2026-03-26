clear;
%close all;

addpath('..');

% topol=[1,2,3;2,3,4;3,4,5];
% coord=[-1,1;0,0;1,1;2,0;3,1];
% geo=geometry();
% geo.init(topol,coord);
% theme(gcf,"light");
% geo.plot(label=[1,2,3]);

% L=1;
% H=2;
% model=createpde();
% rect=[3;4;-L;L;L;-L;-H;-H;H;H];
% geometryFromEdges(model,decsg(rect));
% generateMesh(model,'Hmax',0.1,'GeometricOrder','linear');
% geo=geometry();
% geo.readmodel(model);
% %theme(gcf,"light");
% %geo.plot(label=[1,2,3]);

% R=3;
% C=[0,0];
% model=createpde();
% circ=[1;C(1);C(2);R];
% geometryFromEdges(model,decsg(circ));
% generateMesh(model,'Hmax',1,'GeometricOrder','linear');
% geo=geometry();
% geo.readmodel(model);

h=0.05;

L=1;
H=1;
model=createpde();
rect=[3;4;0;1;1;0;0;0;H;H];
geometryFromEdges(model,decsg(rect));
generateMesh(model,'Hmax',h,'GeometricOrder','linear');
geo=geometry();
geo.readmodel(model);

Omega=geo.whole;
GammaD=Omega.bound;
GammaN=Omega.bound;
Gamma=Omega.skel;
U=space(Omega,2);
R0=space(GammaD,0);
R1=space(GammaN,0);
M=space(Gamma,0);
H=space(GammaN,0);
u=field(U);
r0=field(R0);
r1=field(R1);
m=field(M);

f=field(U);
f.map(@(x) ones(size(x)));
f.eval;

g0=field(R0);
g0.map(@(x) zeros(size(x)));
g0.eval;

g1=field(R1);
g1.map(@(x) zeros(size(x)));
g1.eval;

h0=field(H);
h0.map(@(x) zeros(size(x)));
h0.eval;

h1=field(H);
h1.map(@(x) zeros(size(x)));
h1.eval;

alpha=1;
beta=0;

int=alpha*squeeze(sum(reshape(U.shape{2,3},2,2,U.ref.ndof,1,[]).*reshape(U.shape{2,3},2,2,1,U.ref.ndof,[]),[1,2]));
int=int+beta*squeeze(sum(reshape(U.shape{2,2},2,U.ref.ndof,1,[]).*reshape(U.shape{2,2},2,1,U.ref.ndof,[]),1));
A=asmbmat(int,2,2,2,2,U,U);
l=asmbvec(f.val{2,1}.*U.shape{2,1},2,2,2,U);
int=h0.val{2,1}.*U.shape{1,1};
shapen=squeeze(sum(U.shape{1,2}.*reshape(geo.gauss{2}.trace{1}.normal,2,1,[]),1));
int=int+h1.val{2,1}.*shapen;
l=l+asmbvec(int,2,1,1,U);

int=reshape(U.shape{1,2},2,U.ref.ndof,[]).*reshape(geo.gauss{2}.trace{1}.normal,2,1,[]);
int=squeeze(sum(int,1));
int=reshape(int,U.ref.ndof,1,[]).*reshape(M.shape{2,1},1,M.ref.ndof,[]);
C=asmbmat(int,2,1,1,2,U,M);

int=reshape(U.shape{1,1},U.ref.ndof,1,[]).*reshape(R0.shape{2,1},1,R0.ref.ndof,[]);
B0=asmbmat(int,2,1,1,2,U,R0);
int=g0.val{1,1}.*R0.shape{1,1};
n0=asmbvec(int,1,1,1,R0);

int=reshape(shapen,U.ref.ndof,1,[]).*reshape(R1.shape{2,1},1,R1.ref.ndof,[]);
B1=asmbmat(int,2,1,1,2,U,R1);
int=g1.val{1,1}.*R1.shape{1,1};
n1=asmbvec(int,1,1,1,R1);

mlt=multi;
mlt.add(U);
mlt.add(R0);
mlt.add(R1);
mlt.add(M);
mlt.resblk{1}=l;
mlt.resblk{2}=n0;
mlt.resblk{3}=n1;
mlt.stiffblk{1,1}=A;
mlt.stiffblk{1,2}=B0;
mlt.stiffblk{2,1}=B0';
mlt.stiffblk{1,3}=B1;
mlt.stiffblk{3,1}=B1';
mlt.stiffblk{1,4}=-C;
mlt.stiffblk{4,1}=-C';

stiff=mlt.stiff;
res=mlt.res;

dof=stiff\res;

u.dof=mlt.comp(dof,1);
r0.dof=mlt.comp(dof,2);
r1.dof=mlt.comp(dof,3);
m.dof=mlt.comp(dof,4);

figure;
theme(gcf,"light");
u.eval;
u.plot(trisurf=1,LineStyle="-");
colorbar;