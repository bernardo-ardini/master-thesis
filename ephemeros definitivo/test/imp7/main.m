clear;
clc;
close all;
addpath('../..');
addpath('..');

% model

md=modeloldold();

% geometry

% L=2;
% H=1;
% model=createpde();
% rect=[3;4;-L;L;L;-L;-H;-H;H;H];
% geometryFromEdges(model,decsg(rect));
% h=0.1; % volendo anche h=0.1 funziona, mentre con h più piccolo non più
% generateMesh(model,'Hmax',h,'GeometricOrder','linear');
% geo=geometry();
% geo.readmodel(model);

geo=geometry();
geo.readgmsh("aaa.geo");

figure;
theme(gcf,"light");
geo.plot();
drawnow;


md.geo=geo;
md.name="imp2";

% domains

Omega=geo.whole;
GammaD=Omega.bound.sub(@(x) 0*(abs(abs(x(2,:))-1)<geo.tol));
GammaDn=Omega.bound.sub(@(x) (abs(x(2,:)+1)<geo.tol)|(abs(x(1,:)+10)<geo.tol)|(abs(x(1,:)-10)<geo.tol));

md.Omega=Omega;
md.GammaD=GammaD;
md.GammaDn=GammaDn;

% time

md.T=60;
md.dt0=1.5e-1;
md.dtmax=5;
md.told=Inf;
md.toln=1e-5;
md.tols=1e-4;
md.tolal=Inf;
md.tolpi=Inf;
md.scale=10;

% data

%pl.dat.f=@(x,t) repmat([-1*(1-cos(2*pi*t/10));2*(1-cos(2*pi*t/10))],[1,size(x,2)]);
%pl.dat.f=@(x,t) [(10*((x(2,:)-0.7).*(x(2,:)-1.3)<0)-10*(x(2,:)>1.4))*(1-cos(2*pi*t/10));-0*(1-cos(2*pi*t/10)).*ones(1,size(x,2))];
%pl.dat.f=@(x,t) [10*(x(2,:)>0.5)*(1-cos(2*pi*t/10));-0*(1-cos(2*pi*t/10)).*ones(1,size(x,2))];
md.dat.f=@(x,t) 0*[-0.1*(1-cos(t/pi/2));-0.1*(1-cos(pi*t/2))]*ones(1,size(x,2));
md.dat.g=@(x,t) -4*0.5*(1-cos(pi*t/md.T))*[((abs(x(1,:)-10)<geo.tol));0*x(2,:)];

% material

I=tensorprod(eye(2),eye(2));
Id=reshape(eye(2^2),[2,2,2,2]);
Tr=permute(Id,[1,4,2,3]);
md.mat.C=100*I+100*(Id+Tr);
eta=5;
md.mat.D=eta*I+eta*(Id+Tr);

md.mat.a=4;
md.mat.b=4;
md.mat.c=4;
md.mat.d=eta;
md.mat.e=eta;

c=0.5;
phi=30;
psi=30;

a=sqrt(2)*sind(phi);
b=sqrt(2)*sind(psi);
f0=sqrt(2)*c*cosd(phi);

function [YY,dY]=Y(p,a,b)
    mask=(p>0);
    YY=mask.*(a-b).*p;
    dY=(a-b)*mask.*ones(size(p));
end
md.mat.Y=@(p) Y(p,a,b);

md.mat.R=@(rho) R(rho,1e-3);
md.mat.Q=@(rho,be) Q(rho,be,1e-4);

function [dWd,d2Wd]=Wd(al,f0)
    mask=al<0;
    C=1e6;
    dWd=f0*exp(-al)+C*mask.*al;
    d2Wd=-f0*exp(-al)+C*mask;
end
md.mat.Wd=@(al) Wd(al,f0);

function [LLL,dL,d2L]=LL(al,b)
    LLL=b*(1-exp(-al));
    dL=b*exp(-al);
    d2L=-b*exp(-al);
end
md.mat.L=@(al) LL(al,b);

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

md.init;
%pl.unk.dotz.dof=rand(pl.spc.Z.ndof,1);
%[res,stiff]=pl.assemble;

md.solve;
% 
% field.val=pl.unk.U.hist.proc.val{end};
% field.max=pl.unk.U.hist.proc.max;
% field.min=pl.unk.U.hist.proc.min;
% 
% figure;
% theme(gcf,"light");
% geo.plot("field",field);
% colorbar;