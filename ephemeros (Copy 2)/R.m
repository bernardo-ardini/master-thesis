function [dRR,d2RR]=R(rho,epsilon)
d=size(rho,1);
N=size(rho,3);
r=squeeze(pagenorm(rho,"fro"));

r=r(:);
I=reshape(eye(d^2),[d,d,d,d]);
I=tensorprod(I,ones(N,1));
rhorho=reshape(rho,1,1,d,d,N).*reshape(rho,d,d,1,1,N);

mask1=(r<=epsilon);

dRR=zeros(d,d,N);
d2RR=zeros(d,d,d,d,N);

dRR(:,:,mask1)=rho(:,:,mask1)/epsilon;
dRR(:,:,~mask1)=rho(:,:,~mask1)./reshape(r(~mask1),1,1,[]);

d2RR(:,:,:,:,mask1)=I(:,:,:,:,mask1)/epsilon;
d2RR(:,:,:,:,~mask1)=1./reshape(r(~mask1),1,1,1,1,[]).*(I(:,:,:,:,~mask1)-rhorho(:,:,:,:,~mask1)./reshape(r(~mask1),1,1,1,1,[]).^2);
end

% function [dRR,d2RR]=R(rho,epsilon)
% d=size(rho,1);
% r=pagenorm(rho,"fro");
% mask=(r<=epsilon);
% 
% r=reshape(r,1,1,size(rho,3));
% dRR=zeros(d,d,size(rho,3));
% dRR(:,:,mask)=1/(2*epsilon).*rho(:,:,mask);
% dRR(:,:,~mask)=rho(:,:,~mask)./r(:,:,~mask)-epsilon./(2*r(:,:,~mask).^2).*rho(:,:,~mask);
% 
% r=reshape(r,1,1,1,1,size(rho,3));
% I=reshape(eye(d^2),[d,d,d,d]);
% I=tensorprod(I,ones(size(rho,3),1));
% rhorho=reshape(rho,1,1,d,d,size(rho,3)).*reshape(rho,d,d,1,1,size(rho,3));
% d2RR(:,:,:,:,mask)=I(:,:,:,:,mask)/(2*epsilon);
% d2RR(:,:,:,:,~mask)=((1./r(:,:,:,:,~mask)-epsilon./(2*r(:,:,:,:,~mask).^2)).*I(:,:,:,:,~mask)+(epsilon./r(:,:,:,:,~mask).^4-1./r(:,:,:,:,~mask).^3).*rhorho(:,:,:,:,~mask));
% end

% function [dRR,d2RR]=R(rho)
%     epsilon=1e-3;
% 
%     d=size(rho,1);
%     N=size(rho,3);
% 
%     r=squeeze(pagenorm(rho,"fro"));
%     r=r(:);
%     R=sqrt(r.^2+epsilon^2);
% 
%     I=reshape(eye(d^2),[d,d,d,d]);
%     I=tensorprod(I,ones(size(rho,3),1));
%     rhorho=reshape(rho,1,1,d,d,size(rho,3)).*reshape(rho,d,d,1,1,size(rho,3));
% 
%     dRR=rho./reshape(R,1,1,N);
%     d2RR=reshape(1./R,1,1,1,1,N).*(I-rhorho./reshape(R.^2,1,1,1,1,N));
% end
% md.mat.R=@R;