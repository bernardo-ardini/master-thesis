clear;
close all;

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
% generateMesh(model,'Hmax',1,'GeometricOrder','linear');
% geo=geometry();
% geo.readmodel(model);
% %theme(gcf,"light");
% %geo.plot(label=[1,2,3]);

R=3;
C=[0,0];
model=createpde();
circ=[1;C(1);C(2);R];
geometryFromEdges(model,decsg(circ));
generateMesh(model,'Hmax',0.8,'GeometricOrder','linear');
geo=geometry;
geo.readmodel(model);

Omega=geo.whole;
GammaD=Omega.bound.sub(@(x) x(2,:)>0);
U=space(Omega,-1);
u=field(U);

u.map(@(x) x(1,:));
u.eval;
u.plot(trisurf=1);