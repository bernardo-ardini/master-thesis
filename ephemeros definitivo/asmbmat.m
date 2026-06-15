function out=asmbmat(int,k,j,roi,coi,rospace,cospace)
    int=reshape(int,rospace.ref.ndof,cospace.ref.ndof,[]);
    [H,K,G]=size(int);
    [idh,idk,idg]=ndgrid(1:H,1:K,1:G);
    V=int.*reshape(rospace.geo.gauss{k}.trace{j}.measure,1,1,G);

    R=rospace.ind{roi};
    C=cospace.ind{coi};

    row=R(sub2ind(size(R),idh(:),idg(:)));
    col=C(sub2ind(size(C),idk(:),idg(:)));
    val=V(:);

    row=row(:);
    col=col(:);
    val=val(:);

    mask=(row>0)&(col>0);

    out=sparse(row(mask),col(mask),val(mask),rospace.ndof,cospace.ndof);

    % [H,K,G]=size(int);
    % 
    % R=rospace.ind{roi};
    % C=cospace.ind{coi};
    % 
    % row=zeros(G,H,K);
    % col=zeros(G,H,K);
    % val=zeros(G,H,K);
    % 
    % for g=1:G
    %     [co,ro]=meshgrid(C(:,g),R(:,g));
    %     row(g,:,:)=ro;
    %     col(g,:,:)=co;
    %     val(g,:,:)=int(:,:,g)*rospace.geo.gauss{k}.trace{j}.measure(g);
    % end
    % 
    % row=row(:);
    % col=col(:);
    % val=val(:);
    % 
    % mask=(row>0)&(col>0);
    % 
    % out=sparse(row(mask),col(mask),val(mask),rospace.ndof,cospace.ndof);
end