clear;
close all;

syms E1 E2 G J;

A=1000;
B=10000;
K=0;
l=0.1;
p=2;
alpha=1;
r=2;
assert((p>1)&&(alpha+1<=p));
W=A*((E1-1)^p+(E2-1)^p+l^2*G^p)/J^alpha+B*(J-1)^r+K*(-log(J)+J-1);

J=abs(E1*E2);

W1=diff(subs(W),E1);
W2=diff(subs(W),E2);

W11=diff(subs(W),E1,E1);
W12=diff(subs(W),E1,E2);
W22=diff(subs(W),E2,E2);
WGG=diff(subs(W),G,G);

syms dchi1 chi2 P dP;

C=5;
D=10;
d=50;
q=1.1;
E=0.01;
H=C*(P-1)^2+D*dP^2;
%H=C*(P^q+q/d*exp(-d*(P-1)))+D*dP^2;

HP=diff(H,P);

Y=20;

G=0;
dP=0;

N=100;
tab=zeros(N,10);
ss=zeros(N,5);
inx=zeros(N,1);
dchi1min=1+1e-3;
dchi1max=1.3;

for i=1:N
    lin=linspace(dchi1min,dchi1max,N);
    dchi1=lin(i);

    E1=dchi1;
    E2=chi2;
    han=matlabFunction(subs(subs(W2)));
    tmp=fsolve(han,0.9);
    %tmp=vpasolve(subs(subs(W2==0)),chi2,[0,Inf]);
    solel.chi2=tmp;
    solel.P=1;
    solel.mu=subs(subs(W1*dchi1/P^2-W2*chi2),[chi2,P],[tmp,1]);
    solel.F=subs(subs(W1),[chi2,P],[tmp,1]);
    solel.lambda=subs(subs((-W22+W12^2/W11)/WGG*P^3),[chi2,P],[tmp,1]);

    E1=dchi1/P;
    E2=chi2*P;

    Op=W+H+Y*log(P);
    Om=W+H-Y*log(P);

    ss(i,1)=subs(diff(subs(subs(Om)),P,P),[chi2,P],[tmp,1]);
    ss(i,2)=subs(diff(subs(subs(Op)),P,P),[chi2,P],[tmp,1]);

    han=matlabFunction([subs(W2),subs(W1*dchi1/P^2-W2*chi2-HP-Y/P)],'Vars',{[chi2,P]});
    sol=fsolve(han,[0.9,1.1]);
    solpl.chi2=sol(1);
    solpl.P=sol(2);
    %solpl=vpasolve(subs([subs(W2==0),subs(W1*dchi1/P^2-W2*chi2-HP-Y/P==0)]),[chi2,P],[0,1;0,5]);
    solpl.F=subs(subs(W1)/P,[chi2,P],[solpl.chi2,solpl.P]);
    solpl.lambda=subs(subs((-W22+W12^2/W11)/WGG*P^3),[chi2,P],[solpl.chi2,solpl.P]);

    Op=W+H+Y*log(P/solpl.P);
    Om=W+H-Y*log(P/solpl.P);

    ss(i,3)=subs(diff(subs(subs(Om)),P,P),[chi2,P],[solpl.chi2,solpl.P]);
    ss(i,4)=subs(diff(subs(subs(Op)),P,P),[chi2,P],[solpl.chi2,solpl.P]);

    inx(i)=solel.mu<=Y;
    tab(i,:)=[dchi1,solel.chi2,solel.P,solel.mu,solpl.chi2,solpl.P,solel.F,solpl.F,solel.lambda,solpl.lambda];

    syms dchi1;
    E1=dchi1/P;
    E2=chi2*P;
    he=double(subs(hessian(subs(subs(Op)),[dchi1,chi2,P]),[dchi1,chi2,P],[lin(i),solel.chi2,1]));
    ss(i,5)=min(eig(he));
    he=double(subs(hessian(subs(subs(Op)),[dchi1,chi2,P]),[dchi1,chi2,P],[lin(i),solpl.chi2,solpl.P]));
    [ss(i,6),j]=min(eig(he));

    % if ss(i,6)<=0
    %     [V,D]=eig(he);
    %     disp(V(:,j));
    %     assert(1<0);
    % end
end

figure(3);
%plot(tab(:,1),ss(:,1),"k");
hold on;
%plot(tab(:,1),ss(:,3),"k");
hold on;
%plot(tab(:,1),ss(:,2),"r");
hold on;
plot(tab(:,1),ss(:,4),"r");

figure(4);
%plot(tab(:,1),ss(:,5),"k");
hold on;
plot(tab(:,1),ss(:,6),"k");
hold on;
plot(tab(:,1),0*ss(:,6),"--k");

tmp=tab;

el=find(inx);
pl=find(~inx);

tab=tmp(:,1);
tab(:,2)=[tmp(el,2);tmp(pl,5)];
tab(:,3)=[tmp(el,3);tmp(pl,6)];
tab(:,4)=[tmp(el,7);tmp(pl,8)];
tab(:,5)=[tmp(el,9);tmp(pl,10)];

figure(1);
subplot(2,2,1);
plot(tab(:,1),tab(:,2),"k");
xlabel("dchi1");
ylabel("chi2");
xlim([dchi1min,dchi1max]);
pbaspect([1,1,1]);

subplot(2,2,2);
plot(tab(:,1),tab(:,3),"k");
xlabel("dchi1");
ylabel("P");
xlim([dchi1min,dchi1max]);
pbaspect([1,1,1]);

subplot(2,2,4);
plot(tab(:,1),tab(:,4),"k");
xlabel("dchi1");
ylabel("F");
xlim([dchi1min,dchi1max]);
%ylim([1,30])
pbaspect([1,1,1]);

figure(2);
plot(tab(:,1),tab(:,5),"k");
xlabel("dchi1");
ylabel("\lambda");
xlim([dchi1min,dchi1max]);
pbaspect([1,1,1]);

writematrix(tab(:,[1,4]),"instabilitypl.dat","Delimiter","tab");