A=Q(1:3,1:3);
B=Q(4:5,1:3);
C=Q(4:5,4:5);
q=C-B*(A\B');

ei=eig(A);
disp(ei(1));
disp(ei(2));
disp(ei(3));

[C4,C5]=meshgrid(linspace(-2,2,1000),linspace(-2,2,1000));

Z=q(1,1)*C4.^2+(q(1,2)+q(2,1))*C4.*C5+q(2,2)*C5.^2;

figure;
colormap([0 0 0;1 1 1]); 

imagesc([-2,2],[-2,2],Z>0);