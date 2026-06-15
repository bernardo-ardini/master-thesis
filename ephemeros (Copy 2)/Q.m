function [drhoQ,dbeQ,drhorhoQ,drhobeQ,dbebeQ]=Q(rho,be,epsilon)
d=size(rho,1);
N=size(rho,3);
r=squeeze(pagenorm(rho,"fro"));

be=be(:);
r=r(:);

mask0=(r>1e-12);
%mask1=(be>=r);
mask2=(be(:)<r)&(be(:)>-r);
mask3=(be(:)<=-r);

I=reshape(eye(d^2),[d,d,d,d]);
I=tensorprod(I,ones(N,1));

drhoQ=zeros(d,d,N);
dbeQ=zeros(1,N);
drhorhoQ=zeros(d,d,d,d,N);
drhobeQ=zeros(d,d,N);
dbebeQ=zeros(1,N);

drhoQ(:,:,mask0&mask2)=1/(2*epsilon)*reshape(1-be(mask0&mask2)./r(mask0&mask2),1,1,[]).*rho(:,:,mask0&mask2);
drhoQ(:,:,mask0&mask3)=1/epsilon*rho(:,:,mask0&mask3);

dbeQ(1,(~mask0&(be<=0))|(mask0&mask3))=1/epsilon*be((~mask0&(be<=0))|(mask0&mask3));
dbeQ(1,mask0&mask2)=1/(2*epsilon)*(be(mask0&mask2)-r(mask0&mask2));

drhorhoQ(:,:,:,:,(~mask0&(be<=0))|(mask0&mask3))=1/epsilon*I(:,:,:,:,(~mask0&(be<=0))|(mask0&mask3));
drhorhoQ(:,:,:,:,mask0&mask2)=1/(2*epsilon)*(reshape(1-be(mask0&mask2)./r(mask0&mask2),1,1,1,1,[]).*I(:,:,:,:,mask0&mask2)+reshape(be(mask0&mask2)./r(mask0&mask2).^3,1,1,1,1,[]).*reshape(rho(:,:,mask0&mask2),d,d,1,1,[]).*reshape(rho(:,:,mask0&mask2),1,1,d,d,[]));

drhobeQ(:,:,mask0&mask2)=-1/(2*epsilon)*rho(:,:,mask0&mask2)./reshape(r(mask0&mask2),1,1,[]);

dbebeQ(1,(~mask0&(be<=0))|(mask0&mask3))=1/epsilon*ones(1,nnz((~mask0&(be<=0))|(mask0&mask3)));
dbebeQ(1,mask0&mask2)=1/(2*epsilon)*ones(1,nnz(mask0&mask2));
end

% function [drhoQ,dbeQ,drhorhoQ,drhobeQ,dbebeQ]=Q(rho,be)
%     epsilon=13;
%     mu=12;
% 
%     d=size(rho,1);
%     N=size(rho,3);
%     r=squeeze(pagenorm(rho,"fro"));
% 
%     be=be(:);
%     r=r(:);
% 
%     I=reshape(eye(d^2),[d,d,d,d]);
%     I=repmat(I,[1,1,1,1,N]);
%     rhorho=reshape(rho,d,d,1,1,N).*reshape(rho,1,1,d,d,N);
% 
%     R=sqrt(r.^2+epsilon^2);
%     dist=R-be;
% 
%     drhoQ=reshape(-mu./dist./R,1,1,N).*rho;
%     dbeQ=reshape(mu./dist,1,N);
% 
%     drhorhoQ=reshape(mu./dist./R,1,1,1,1,N).*(reshape((1./dist+1./R)./R,1,1,1,1,N).*rhorho-I);
%     drhobeQ=-reshape(mu./dist.^2./R,1,1,[]).*rho;
%     dbebeQ=reshape(mu./dist.^2,1,N);
% end