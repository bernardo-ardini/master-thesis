clear;
close all;
clc;

syms E1 E2 E3 P dP P0 real;

A=1.1e6;B=1.1e6;l=4;p=2;r=2;alpha=1;
assert((alpha>0)&(r>1)&(alpha<=p-1)&(p>1+alpha/r));
W=A*(((E1-1)^2)^(p/2)+((E2-1)^2)^(p/2)+((l*E3)^2)^(p/2))/(E1*E2)^alpha+B*((E1*E2-1)^2)^(r/2);

C=0.3e3;D=10;E=1e3;l=0.1;Pi=0.1;q=1.3;delta=17;
%H=C*((P-1)^2)^(q/2)+D/(P-Pi)+E*((l*dP)^2)^(q/2);
H=C*(P^q+q/delta*exp(-delta*(P-1)))+D/(P-Pi)+E*(l*dP)^2;

Y=3.8e3;
ep=1e-9;
%D=Y*sqrt(ep+log(P/P0)^2);
D=Y*log(P/P0);

syms eta y1 dy1 y2 dy2 y3 dy3 real;

in=[E1,E2,E3,P,dP,P0];
out=[eta/y3,y2*y3,0,y3,0,y3];

W0=subs(W,in(1:end-1),out(1:end-1));
H0=subs(H,in(1:end-1),out(1:end-1));
D0=subs(D,in(1:end-1),out(1:end-1));
W1=subs(diff(W,E1),in,out);
W2=subs(diff(W,E2),in,out);
W11=subs(diff(W,E1,E1),in,out);
W12=subs(diff(W,E1,E2),in,out);
W22=subs(diff(W,E2,E2),in,out);
W33=subs(diff(W,E3,E3),in,out);
H1=subs(diff(H,P),in,out);
H11=subs(diff(H,P,P),in,out);
H22=subs(diff(H,dP,dP),in,out);
D2=subs(diff(D,P),in,out);
D22=subs(diff(D,P,P),in,out);

in=[E1,E2,E3,P,dP];
out=[eta*dy1/y3,y2*y3,dy2/y3,y3,dy3];

W=subs(W,in,out);
H=subs(H,in,out);
D=subs(D,in,out);

in=[dy1,dy2,dy3,P0];
out=[1,0,0,y3];

A=subs(hessian(W+H+D,[dy1,dy2,dy3]),in,out);
B=subs(jacobian(gradient(W+H+D,[dy1,dy2,dy3]),[y1,y2,y3]),in,out);
C=subs(hessian(W+H+D,[y1,y2,y3]),in,out);

syms x L delta real;

h=l;
V=[delta*(1-exp(-x/delta))-x/L*delta*(1-exp(-L/delta)),0,0;sin(pi*x/L),0,0;0,cos(pi/L*x),0;0,1,0;0,-x/delta*exp(-x/delta),0;0,-exp(-x/delta),0;0,0,exp(-x/delta)]';
dV=diff(V,x);
Q=int(dV'*A*dV+dV'*B*V+V'*B'*dV+V'*C*V,x,[0,L]);

N=20;
res=cell(N,1);
etas=linspace(1,1.003,N);
etas=[etas(1:end-1),linspace(1.003,1.35,N)];
N=length(etas);
yield=0;
L0=6e5;

% load deQQ.mat;
% syms h;

for i=1:N
    if yield==0
        han=matlabFunction(subs(W2,[eta,y3],[etas(i),1]));
        sol=fsolve(han,0.9);
        ma=subs(subs(eta/y3^2*W1-y2*W2-H1),[eta,y3,y2],[etas(i),1,sol]);
        res{i}.y2=sol;
        res{i}.y3=1;
        res{i}.yield=0;
        res{i}.F=double(subs(subs(W1/y3),[eta,y3,y2],[etas(i),1,sol]));
        if ma>Y
            yield=1;
        end
    end
        
    if yield==1
        han=matlabFunction(subs([W2,eta/y3^2*W1-y2*W2-H1-D2],eta,etas(i)),'Vars',{[y2,y3]});
        sol=fsolve(han,[0.9,1.1]);

        in=[eta,y2,y3,L];
        out=[etas(i),sol(1),sol(2),L0];

        q=subs(Q,in,out);
        %res{i}.delta=vpasolve(diff(det(q),delta),delta);
        deltas=[l,10*l,L0/10,L0/3,L0/2,2/3*L0,L0,10*L0];

        for k=1:length(deltas)
            q=double(subs(subs(Q,in,out),delta,deltas(k)));
            [V,D]=eigs(q,1,'smallestreal');
            res{i}.Q{k}=q;
            res{i}.V{k}=V;
            res{i}.D(k)=D;
        end
        %han=matlabFunction(diff(det(q),delta));
        %fsolve(han,l);
        %gradient(det(q),[a];

        % [V,D]=eigs(q,1,'smallestreal');
        % res{i}.V=V;
        % res{i}.D=D;
        % res{i}.Q=q;
    
        res{i}.y2=sol(1);
        res{i}.y3=sol(2);
        res{i}.yield=1;
        res{i}.F=double(subs(subs(W1/y3),[eta,y2,y3],[etas(i),sol(1),sol(2)]));
    end
end

figure;
theme(gcf,"light");
plot(etas,arrayfun(@(i) res{i}.y2,1:N));
title("tranversal deformation");

figure;
theme(gcf,"light");
plot(etas,arrayfun(@(i) res{i}.y3,1:N));
title("plastic distortion");

figure;
theme(gcf,"light");
plot(etas,arrayfun(@(i) res{i}.F,1:N));
pbaspect([1 1 1])
title("force");

ind=1:N;
ind=ind(arrayfun(@(i) res{i}.yield,1:N)==1);

% figure;
% theme(gcf,"light");
% theme(gcf,"light");
% plot(etas(ind),arrayfun(@(i) res{i}.delta,ind));
% pbaspect([1 1 1])
% title("delta");

figure;
theme(gcf,"light");
etan=zeros(1,length(deltas));
for k=1:length(deltas)
    plot(etas(ind),arrayfun(@(i) res{i}.D(k),ind),"LineWidth",2);
    hold on;
    han=@(t) interp1(etas(ind),arrayfun(@(i) res{i}.D(k),ind),t,'spline');
    etan(k)=fzero(han,1);
end
legend(string(deltas));
pbaspect([1 1 1]);
title("det Q");

figure;
theme(gcf,"light");
plot(etas(ind),arrayfun(@(i) res{i}.V{3}(end),ind),"-","LineWidth",2);
hold on;
plot(etas(ind),arrayfun(@(i) res{i}.V{3}(end-1),ind),"--","LineWidth",2);
% hold on;
% plot(etas(ind),arrayfun(@(i) res{i}.V{3}(3),ind),"LineWidth",2);
legend("1","2","3");
pbaspect([1 1 1]);
title("det Q");

% figure;
% theme(gcf,"light");
% plot(etas(ind),arrayfun(@(i) res{i}.V(4),ind),"b");
% hold on;
% plot(etas(ind),arrayfun(@(i) res{i}.V(5),ind),"r");
% pbaspect([1 1 1])
% title("components of the eigenvector");

han=@(t)-interp1(etas,arrayfun(@(i) res{i}.F,1:N),t,'spline');
[etam,~]=fminbnd(han,1,2);
etay=mean([etas(ind(1)-1),etas(ind(1))]);
% han=@(t) interp1(etas(ind),arrayfun(@(i) res{i}.D,ind),t,'spline');
% etan=fzero(han,1);

% Q=res{end}.Q;
% A=Q(1:3,1:3);
% B=Q(4:5,1:3);
% C=Q(4:5,4:5);
% q=C-B*(A\B');
% 
% [C4,C5]=meshgrid(linspace(-2,2,1000),linspace(-2,2,1000));
% 
% Z=q(1,1)*C4.^2+(q(1,2)+q(2,1))*C4.*C5+q(2,2)*C5.^2;
% 
% figure;
% colormap([0 0 0;1 1 1]); 
% 
% imagesc([-2,2],[-2,2],Z>0);