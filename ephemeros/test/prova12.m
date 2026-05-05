clear;
clc;
close all;
addpath('..');

% model

pl=plasticity();

% geometry

L=3;
H=1;
model=createpde();
rect=[3;4;-L;L;L;-L;-H;-H;H;H];
geometryFromEdges(model,decsg(rect));
h=0.26;
generateMesh(model,'Hmax',h,'GeometricOrder','linear');
geo=geometry();
geo.readmodel(model);

figure;
theme(gcf,"light");
geo.plot();
drawnow;

pl.geo=geo;

% domains

Omega=geo.whole;
GammaD=Omega.bound.sub(@(x) (abs(abs(x(1,:))-L)<geo.tol));
GammaDn=Omega.bound.sub(@(x) 0*((abs(abs(x(1,:))-L)<geo.tol))+(abs(x(2,:)-H)<geo.tol));

pl.Omega=Omega;
pl.GammaD=GammaD;
pl.GammaDn=GammaDn;

% time

pl.T=2;
pl.setting.step.dtmax=0.1;
pl.setting.step.dtmin=5e-5;
pl.setting.step.dt0=0.01;
pl.setting.step.nmin=15;
pl.setting.step.nmax=4;

% newton-raphson setting

pl.setting.newton.tol=1e-7;
pl.setting.newton.maxiter=30;
pl.setting.step.maxgood=4;
pl.setting.tol=0.1;

% data

%pl.dat.f=@(x,t) repmat([-1*(1-cos(2*pi*t/10));2*(1-cos(2*pi*t/10))],[1,size(x,2)]);
%pl.dat.f=@(x,t) [(10*((x(2,:)-0.7).*(x(2,:)-1.3)<0)-10*(x(2,:)>1.4))*(1-cos(2*pi*t/10));-0*(1-cos(2*pi*t/10)).*ones(1,size(x,2))];
%pl.dat.f=@(x,t) [10*(x(2,:)>0.5)*(1-cos(2*pi*t/10));-0*(1-cos(2*pi*t/10)).*ones(1,size(x,2))];
pl.dat.f=@(x,t) [-0.00*(1-cos(pi*t/2));-0.5*(1-cos(pi*t/2))]*ones(1,size(x,2));
pl.dat.g=@(x,t) 0*x;
pl.unit.strain=1e-3;

% material

pl.mat.eta=1;
pl.mat.l=0.1;

pl.mat.D.mu=100;
pl.mat.D.la=100;
pl.mat.A=100;

pl.mat.C.mu=400;
pl.mat.C.la=400;

pl.mat.B=0.1;

% pl.mat.l=0.4;
% 
% pl.mat.D.mu=0.0001;
% pl.mat.D.la=0.0001;
% pl.mat.A=0.01;
% pl.mat.E=0.0001;
% 
% pl.mat.C.mu=400;
% pl.mat.C.la=400;
% 
% pl.mat.B=0.1;
% pl.mat.F=0.0001;

function [Y,dY]=yield(p)
    alpha=8;
    Y0=1;

    % Y=Y0*zeros(size(p));
    % dY=Y0*zeros(size(p));

    mask=(p>-Y0/alpha);
    Y=mask.*(Y0+alpha*p);
    dY=alpha*mask.*ones(size(p));
end
pl.mat.yield=@yield;

function [dG,d2G]=diss(dotz)
    c=1e-4;

    d=size(dotz,1);
    r=pagenorm(dotz,"fro");
    mask=(r<=c);

    r=reshape(r,1,1,size(dotz,3));
    dG=zeros(d,d,size(dotz,3));
    dG(:,:,mask)=1/(2*c).*dotz(:,:,mask);
    dG(:,:,~mask)=dotz(:,:,~mask)./r(:,:,~mask)-c./(2*r(:,:,~mask).^2).*dotz(:,:,~mask);

    r=reshape(r,1,1,1,1,size(dotz,3));
    I=reshape(eye(d^2),[d,d,d,d]);
    I=tensorprod(I,ones(size(dotz,3),1));
    dotzdotz=reshape(dotz,1,1,d,d,size(dotz,3)).*reshape(dotz,d,d,1,1,size(dotz,3));
    d2G(:,:,:,:,mask)=I(:,:,:,:,mask)/(2*c);
    d2G(:,:,:,:,~mask)=((1./r(:,:,:,:,~mask)-c./(2*r(:,:,:,:,~mask).^2)).*I(:,:,:,:,~mask)+(c./r(:,:,:,:,~mask).^4-1./r(:,:,:,:,~mask).^3).*dotzdotz(:,:,:,:,~mask));
end
pl.mat.diss=@diss;

% function [dG,d2G]=diss(dotz)
%     c=1e-6;
% 
%     d=size(dotz,1);
%     r=pagenorm(dotz,"fro");
%     mask=(r<=c);
% 
%     r=reshape(r,1,1,size(dotz,3));
%     dG=zeros(d,d,size(dotz,3));
%     dG(:,:,mask)=1/c.*dotz(:,:,mask);
%     dG(:,:,~mask)=dotz(:,:,~mask)./r(:,:,~mask);
% 
%     r=reshape(r,1,1,1,1,size(dotz,3));
%     I=reshape(eye(d^2),[d,d,d,d]);
%     I=tensorprod(I,ones(size(dotz,3),1));
%     dotzdotz=reshape(dotz,1,1,d,d,size(dotz,3)).*reshape(dotz,d,d,1,1,size(dotz,3));
%     d2G(:,:,:,:,mask)=I(:,:,:,:,mask)/c;
%     d2G(:,:,:,:,~mask)=1./r(:,:,:,:,~mask).*I(:,:,:,:,~mask)-1./r(:,:,:,:,~mask).^3.*dotzdotz(:,:,:,:,~mask);
% end
% pl.mat.diss=@diss;

% init

pl.init;
%pl.unk.dotz.dof=rand(pl.spc.Z.ndof,1);
%[res,stiff]=pl.assemble;

pl.solve;
% 
% field.val=pl.unk.U.hist.proc.val{end};
% field.max=pl.unk.U.hist.proc.max;
% field.min=pl.unk.U.hist.proc.min;
% 
% figure;
% theme(gcf,"light");
% geo.plot("field",field);
% colorbar;