epsilon=1e-6;
pl.setting.step.dt=0.1;
pl.t=0;
%rng(0);

pl.mlt.scale=ones(size(pl.mlt.scale));

dotu0=rand(pl.spc.U.ndof,1);
u0=rand(pl.spc.U.ndof,1);
r0=rand(pl.spc.R.ndof,1);
rn0=rand(pl.spc.Rn.ndof,1);
m0=rand(pl.spc.M.ndof,1);
dotz0=rand(pl.spc.Z.ndof,1);
z0=rand(pl.spc.Z.ndof,1);
la0=rand(pl.spc.L.ndof,1);
p0=rand(pl.spc.P.ndof,1);

pl.unk.dotu.dof=dotu0;
pl.unk.u.dof=u0+pl.dt*pl.unk.dotu.dof;
pl.unk.r.dof=r0;
pl.unk.rn.dof=rn0;
pl.unk.m.dof=m0;
pl.unk.dotz.dof=dotz0;
pl.unk.z.dof=z0+pl.dt*pl.unk.dotz.dof;
pl.unk.la.dof=la0;
pl.unk.p.dof=p0;

[res,stiff]=pl.assemble;

ddotu=1*epsilon*rand(pl.spc.U.ndof,1);
dr=1*epsilon*rand(pl.spc.R.ndof,1);
drn=1*epsilon*rand(pl.spc.Rn.ndof,1);
dm=1*epsilon*rand(pl.spc.M.ndof,1);
ddotz=1*epsilon*rand(pl.spc.Z.ndof,1);
dla=1*epsilon*rand(pl.spc.L.ndof,1);
dp=1*epsilon*rand(pl.spc.P.ndof,1);

delta=[ddotu;dr;dm;ddotz;dla;dp;drn];

pl.unk.dotu.dof=pl.unk.dotu.dof+ddotu;
pl.unk.u.dof=u0+pl.dt*pl.unk.dotu.dof;
pl.unk.r.dof=pl.unk.r.dof+dr;
pl.unk.rn.dof=pl.unk.rn.dof+drn;
pl.unk.m.dof=pl.unk.m.dof+dm;
pl.unk.dotz.dof=pl.unk.dotz.dof+ddotz;
pl.unk.z.dof=z0+pl.dt*pl.unk.dotz.dof;
pl.unk.la.dof=pl.unk.la.dof+dla;
pl.unk.p.dof=pl.unk.p.dof+dp;

[res1,stiff1]=pl.assemble;

dres=res1-res;
dresex=stiff*delta;

format shortE

%disp(norm(dres-dresex,Inf))

for i=1:7
    resi=pl.mlt.comp(res,i);
    res1i=pl.mlt.comp(res1,i);
    dresi=pl.mlt.comp(dres,i);
    dresexi=pl.mlt.comp(dresex,i);
    %disp([resi,res1i,dresi,dresexi,dresexi-dresi])
    disp(norm(dresexi-dresi,Inf));
end

%disp([res,res1,dres,dresex])