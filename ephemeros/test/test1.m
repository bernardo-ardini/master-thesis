clear;
close all;

addpath('..');

topol=[1,2,3;2,3,4;3,4,5];
coord=[-1,1;0,0;1,1;2,0;3,1];
geo=geometry();
geo.init(topol,coord,gord=0);

% L=1;
% H=2;
% model=createpde();
% rect=[3;4;-L;L;L;-L;-H;-H;H;H];
% geometryFromEdges(model,decsg(rect));
% generateMesh(model,'Hmax',0.05,'GeometricOrder','linear');
% geo=geometry();
% geo.readmodel(model);

Omega=geo.whole;
Gamma=Omega.bound;

U=space(Omega,2,type="Vec");
u=field(U);
u.map(@(x) x);
u.eval;

R=space(Gamma,0,type="Vec");
r=field(R);
r.map(@(x) x);
r.eval;