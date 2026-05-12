clear;
close all;
clc;

syms E1 E2 E3 P dP P0 real;

A=1.1e6;B=1.1e6;l=4;p=2;r=2;alpha=1;
assert((alpha>0)&(r>1)&(alpha<=p-1)&(p>1+alpha/r));
W=A*(((E1-1)^2)^(p/2)+((E2-1)^2)^(p/2)+((l*E3)^2)^(p/2))/(E1*E2)^alpha+B*((E1*E2-1)^2)^(r/2);

C=1.3e3;D=10;E=1e3;l=0.1;Pi=0.1;q=1.3;delta=17;
%H=C*((P-1)^2)^(q/2)+D/(P-Pi)+E*((l*dP)^2)^(q/2);
H=C*(P^q+q/delta*exp(-delta*(P-1)))+D/(P-Pi)+E*(l*dP)^2;

Y=3.8e3;
ep=1e-9;
%D=Y*sqrt(ep+log(P/P0)^2);
D=Y*log(P/P0);

syms eta y2 y3 real;

in=[E1,E2,E3,P,dP,P0];
out=[eta/y3,y2*y3,0,y3,0,y3];

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

N=40;
res=cell(N,1);
etas=linspace(1,1.003,N);
etas=[etas(1:end-1),linspace(1.003,1.3,N)];
N=length(etas);
L=201;
%L=30e20;
yield=0;
hs=linspace(0,0.9*L,2);

load deQQ.mat;
syms h;

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

        in=[eta,y2,y3];
        out=[etas(i),sol(1),sol(2)];

        A=double(subs([eta^2*W11/y3^2,0,0;0,W33/y3^2,0;0,0,H22],in,out));
        B=zeros(3,3);
        B(1,2)=double(subs(eta*W12,in,out));
        B(2,1)=-B(1,2);
        B(1,3)=double(subs(eta*y2/y3*W12-eta^2*W11/y3^3-eta*W1/y3^2,in,out));
        B(3,1)=-B(1,3);
        C=zeros(3,3);
        C(2,2)=double(subs(-y3^2*W22,in,out));
        C(2,3)=double(subs(-W2+eta/y3*W12-y2*y3*W22,in,out));
        C(3,2)=C(2,3);
        C(3,3)=double(subs(-eta^2*W11/y3^4-y2^2*W22-2*eta*W1/y3^3+2*eta*y2/y3^2*W12-H11-D22,in,out));

        % syms la;
        % P=det(la^2*A+la*B+C);
        % sop=solve(P==0,la);
        % vp=zeros(3,6);
        % lap=zeros(6,1);
        % for j=1:6
        %     LLL=subs(la^2*A+la*B+C,la,sop(j));
        %     vp(:,j)=null(LLL);
        %     lap(j)=double(sop(j));
        % end
        % 
        % syms la;
        % P=det(la^2*A(1:2,1:2)+la*B(1:2,1:2)+C(1:2,1:2));
        % so0=solve(P==0,la);
        % v0=zeros(2,4);
        % la0=zeros(4,1);
        % for j=1:4
        %     LLL=subs(la^2*A(1:2,1:2)+la*B(1:2,1:2)+C(1:2,1:2),la,so0(j));
        %     v0(:,j)=null(LLL);
        %     la0(j)=double(sop(j));
        % end

        % Mp=[zeros(3,3),eye(3,3);-A\C,-A\B];
        % M0=[zeros(2,2),eye(2,2);-A(1:2,1:2)\C(1:2,1:2),-A(1:2,1:2)\B(1:2,1:2)];
        % 
        % I1p=eye(6,6);
        % I1p(6,:)=[];
        % I2p=eye(6,6);
        % I2p(:,[1,5,6])=[];
        % 
        % I10=[1,0,0,0;0,1,0,0;0,0,0,0;0,0,1,0;0,0,0,1];
        % I20=eye(4,4);
        % I20(:,[1,4])=[];
        % 
        % res{i}.deQ=zeros(length(hs),1);
        % for j=1:length(hs)
        %     [Vp,Dp]=eig(Mp*hs(j));
        %     [V0,D0]=eig(M0*(hs(j)-L));
        %     Q=[I1p*Vp*exp(Dp)/Vp*I2p,-I10*V0*exp(D0)/V0*I20];
        %     res{i}.deQ(j)=det(Q);
        % end

        a1=A(1,1);
        a2=A(2,2);
        a3=A(3,3);
        b1=B(1,2);
        b2=B(1,3);
        c1=C(2,2);
        c2=C(2,3);
        c3=C(3,3);

        syms v10(x) v20(x);
        dv20=diff(v20,x);
        dx=dsolve([a1*diff(v10,x,x)+b1*diff(v20,x)==0,a2*diff(v20,x,x)-b1*diff(v10,x)+c1*v20==0],[v10(L)==0,dv20(L)==0]);

        syms v1p(x) v2p(x) v3p(x);
        dv2p=diff(v2p,x);
        dv3p=diff(v3p,x);
        sx=dsolve([a1*diff(v1p,x,x)+b1*diff(v2p,x)+b2*diff(v3p,x)==0,a2*diff(v2p,x,x)-b1*diff(v1p,x)+c1*v2p+c2*v3p==0,a3*diff(v3p,x,x)-b2*diff(v1p,x)+c2*v2p+c3*v3p==0],[v1p(0)==0,dv2p(0)==0,dv3p(0)==0]);

        syms C1 C2 C3;
        res{i}.deQ=zeros(length(hs),1);

        for j=1:length(hs)
            % matdx=equationsToMatrix(subs([dx.v10;dx.v20;0;diff(dx.v10,x);diff(dx.v20,x)],x,hs(j)),[C1,C2]);
            % matsx=equationsToMatrix(subs([sx.v1p;sx.v2p;sx.v3p;diff(sx.v1p,x);diff(sx.v2p,x)],x,hs(j)),[C1,C2,C3]);
            % Q=[matdx,-matsx];
            res{i}.deQ(j)=double(subs(subs(de,h,hs(j))));
        end

        K=4;
        res{i}.de=zeros(K+1,1);

        res{i}.de(1)=det([c1,c2;c2,c3]);

        for k=1:K
            gamma=pi/L*k;
            M=zeros(3,3);
            M(1,1)=-a1*gamma^2;
            M(1,2)=-b1*gamma;
            M(1,3)=-b2*gamma;
            M(2,1)=M(1,2);
            M(2,2)=-a2*gamma^2+c1;
            M(2,3)=c2;
            M(3,1)=M(1,3);
            M(3,2)=c2;
            M(3,3)=-a3*gamma^2+c3;
            res{i}.de(k+1)=det(M);
        end
    
        res{i}.y2=sol(1);
        res{i}.y3=sol(2);
        res{i}.yield=1;
        res{i}.f=a1*(c1*c3-c2^2)+b1^2*c3+b2^2*c1-2*b1*b2*c2;
        res{i}.F=double(subs(subs(W1/y3),[eta,y2,y3],[etas(i),sol(1),sol(2)]));
    end
end

figure;
plot(etas,arrayfun(@(i) res{i}.y2,1:N));

figure;
plot(etas,arrayfun(@(i) res{i}.y3,1:N));

figure;
plot(etas,arrayfun(@(i) res{i}.F,1:N));
pbaspect([1 1 1])

ind=1:N;
ind=ind(arrayfun(@(i) res{i}.yield,1:N)==1);

han=@(t)-interp1(etas(ind),arrayfun(@(i) res{i}.f,ind),t,'spline');
etaf=fzero(han,1);

han=@(t)-interp1(etas,arrayfun(@(i) res{i}.F,1:N),t,'spline');
[etam,~]=fminbnd(han,1,2);

etay=mean([etas(ind(1)-1),etas(ind(1))]);
etan=zeros(K+1,1);

for k=1:K+1
    figure;
    plot(etas(ind),arrayfun(@(i) res{i}.de(k),ind));
    hold on;

    han=@(t) interp1(etas(ind),arrayfun(@(i) res{i}.de(k),ind),t,'spline');
    etan(k)=fzero(han,1);
end

figure;
plot(etas(ind),arrayfun(@(i) res{i}.f,ind))

figure(3);