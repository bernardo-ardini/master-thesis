clear;
close all;

L=5;
R=0.5;

N=10;
X=linspace(0,L,N)';

eta=1;
zeta=1/eta;

K=0;

x=eta*X+K*L/pi*sin(pi*X/L);
y=R*(zeta-K*cos(pi*X/L));

writematrix([x,y],"instabilityelref.dat","Delimiter","tab")