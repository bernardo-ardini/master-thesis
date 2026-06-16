clear;
close all;
clc;

d=3;

syms E1 E2 E3 P1 dP1 P10 P2 dP2 P20 real;

A=1.1e6;B=1.1e6;l=4;p=2;r=2;alpha=1;
%syms A B l p r alpha real;
%assert((alpha>0)&(r>1)&(alpha<=p-1)&(p>1+alpha/r));
W=A*(((E1-1)^2)^(p/2)+((E2-1)^2)^(p/2)+((l*E3)^2)^(p/2))/(E1*E2^(d-1))^alpha+B*((E1*E2^(d-1)-1)^2)^(r/2);

C=0.1e3;D=10;E=1e3;l=0.01;Pi=0.1;q=1.3;delta=17;
%syms C D E ll Pi q delta real;
%H=C*((P-1)^2)^(q/2)+D/(P-Pi)+E*((l*dP)^2)^(q/2);
H=C*(P^q+q/delta*exp(-delta*(P-1)))+D/(P-Pi)+E*(l*dP)^2;
%H=C*exp()

Y=3.8e3;
%syms Y;
%D=Y*sqrt(ep+log(P/P0)^2);
D=Y*log(P/P0);

syms eta y1 dy1 y2 dy2 y3 dy3 y4 dy4 real;

in=[E1,E2,E3,P1,dP1,P10,P2,dP2,P20];
out=[eta/y3,y2*y3^(1/(d-1)),-y2*y3^((2-d)/(d-1))*y4,y3,0,y3,y4,0,y4];

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
out=[eta*dy1/y3,y2*y3^(1/(d-1)),dy2/y3,y3,dy3];

W=subs(W,in,out);
H=subs(H,in,out);
D=subs(D,in,out);

in=[dy1,dy2,dy3,P0];
out=[1,0,0,y3];

A=subs(hessian(W+H+D,[dy1,dy2,dy3]),in,out);
B=subs(jacobian(gradient(W+H+D,[dy1,dy2,dy3]),[y1,y2,y3]),in,out);
C=subs(hessian(W+H+D,[y1,y2,y3]),in,out);

syms x L real;

V=[sin(pi*x/L),0,0;0,1,0;0,cos(pi*x/L),0;0,0,1;0,0,cos(pi*x/L)]';
dV=diff(V,x);
Q=int(dV'*A*dV+dV'*B*V+V'*B'*dV+V'*C*V,x,[0,L]);

N=40;
res=cell(N,1);
etas=linspace(1,1.003,N);
etas=[etas(1:end-1),linspace(1.003,1.35,N)];
N=length(etas);
yield=0;

% load deQQ.mat;
% syms h;

for i=1:N
    if yield==0
        han=matlabFunction(subs(W2,[eta,y3],[etas(i),1]));
        sol=fsolve(han,0.9);
        ma=subs(subs(eta/y3^2*W1-1/(d-1)*y2*y3^(1/(d-1)-1)*W2-H1),[eta,y3,y2],[etas(i),1,sol]);
        res{i}.y2=sol;
        res{i}.y3=1;
        res{i}.yield=0;
        res{i}.F=double(subs(subs(W1/y3),[eta,y3,y2],[etas(i),1,sol]));
        if ma>Y
            yield=1;
        end
    end
        
    if yield==1
        han=matlabFunction(subs([W2,eta/y3^2*W1-1/(d-1)*y2*y3^(1/(d-1)-1)*W2-H1-D2],eta,etas(i)),'Vars',{[y2,y3]});
        sol=fsolve(han,[0.9,1.1]);

        in=[eta,y2,y3,L];
        out=[etas(i),sol(1),sol(2),6e4];

        q=double(subs(Q,in,out));
        [V,D]=eigs(q,1,'smallestreal');
        res{i}.V=V;
        res{i}.D=D;
        res{i}.Q=q;
    
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

figure;
theme(gcf,"light");
theme(gcf,"light");
plot(etas(ind),arrayfun(@(i) res{i}.D,ind));
pbaspect([1 1 1])
title("smallest eigenvalue");

figure;
theme(gcf,"light");
plot(etas(ind),arrayfun(@(i) res{i}.V(4),ind),"b");
hold on;
plot(etas(ind),arrayfun(@(i) res{i}.V(5),ind),"r");
pbaspect([1 1 1])
title("components of the eigenvector");

han=@(t)-interp1(etas,arrayfun(@(i) res{i}.F,1:N),t,'spline');
[etam,~]=fminbnd(han,1,2);
etay=mean([etas(ind(1)-1),etas(ind(1))]);
han=@(t) interp1(etas(ind),arrayfun(@(i) res{i}.D,ind),t,'spline');
etan=fzero(han,1);

Q=res{60}.Q;
A=Q(1:3,1:3);
B=Q(4:5,1:3);
C=Q(4:5,4:5);
q=C-B*(A\B');

disp(eig(A));

[C4,C5]=meshgrid(linspace(-2,2,1000),linspace(-2,2,1000));

Z=q(1,1)*C4.^2+(q(1,2)+q(2,1))*C4.*C5+q(2,2)*C5.^2;

figure;
colormap([0 0 0;1 1 1]); 

imagesc([-2,2],[-2,2],Z>0);