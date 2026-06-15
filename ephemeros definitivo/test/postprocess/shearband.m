t=34;

[~,s]=min(abs(times-t));
h=md.hist{s};

ccoord=geo.gauss{2}.trace{2}.coord;
vals=h.al{2,1};
vals=vals(:);
x=ccoord(:,1);
y=ccoord(:,2);

x0=[0,-1];
x1=[2,1];
interp=linspace(0,1,1e3)';
gcoord=x0.*(1-interp)+x1.*interp;

al=scatteredInterpolant(x,y,vals,'linear','linear');
vals=al(gcoord(:,1),gcoord(:,2));
vals=smoothdata(vals);

figure;
plot(norm(x1-x0)*interp,vals);