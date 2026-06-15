classdef space < handle
    properties
        geo % geometry
        dom % domain
        p % abs(p) is the poly degree, p>0 continuous, p<1 discontinuous across elements
        m % dim of the element (es. m=2 triangle)
        coord % coord(i,a) coordinate i of node a
        topol % topol(e,a) ath node of element e
        hdshape % shape{r}(I1,...,IN,i1,...,ir,e,g,h) is rth derivative of the component (I1,...,IN) of the hth shape function (h is the index of the dofs) in element e of dim m-1 funct along dirs i1,...,ik computed at local gauss point g
        shape % whshape{k,r}(I1,...,IN,i1,...,ir,h,G)
        num % number of nodes
        ndof % number of dofs
        hdind % hdind(e,h) gives the global dofs index given the local dofs index h in element e
        ind % ind(h,G)
        ref % ref.coord coordinates of nodes in reference element, ref.shape{k}(i1,...,ik,g,a is kth derivative of ath shape funct in reference element along dirs i1,...,ik computed at local gauss point g, ref.num number of nodes in the reference element, ref.ndof number of local dofs
        gauss % gauss.num number of gauss point in reference element, gauss.coord coordinates of gauss points in reference element
        basis % basis.dim is dim N of the type space, basis.B(I1,...,IN,q) is the component (I1,...,IN) of the qth basis vector of the type space
    end
    methods
        function enum(obj)
            obj.num=size(obj.coord,1);
            obj.ref.num=size(obj.ref.coord,1);
            obj.ref.ndof=obj.ref.num*obj.basis.dim;
            if obj.p<=0
                obj.ndof=obj.dom.num*obj.ref.ndof;
            elseif obj.p>=0
                obj.ndof=obj.num*obj.basis.dim;
            end
        end

        function eind(obj)
            if obj.p<=0
                id=1:obj.ndof;
                obj.hdind=reshape(id,obj.dom.num,obj.ref.ndof);
            elseif obj.p>=1
                id=1:obj.ndof;
                id=reshape(id,[obj.num,obj.basis.dim]);
                id=id(obj.topol(:),:);
                id=reshape(id,[obj.dom.num,obj.ref.ndof]);
                obj.hdind=id;
            end

            obj.ind=cell(obj.geo.m,1);
            id=zeros(obj.geo.num(obj.m+1),obj.ref.ndof);
            id(obj.dom.ind,:)=obj.hdind;
            id=id';

            for j=1:obj.geo.m
                if j<=obj.m
                    obj.ind{j}=id(:,obj.geo.gauss{obj.m}.trace{j}.elem(:,1));
                else
                    obj.ind{j}=id(:,obj.geo.gauss{j}.trace{obj.m}.elem(:,2));
                end
            end
        end

        function enode(obj)
            if obj.p==0
                n=1/(obj.m+1)*ones(obj.m+1,1);
            else
                a=zeros(obj.m+1,1);
                a(1)=1;
                A=[-ones(1,obj.m);eye(obj.m)];
                n=a+A*obj.ref.coord'/abs(obj.p);
            end

            Z=reshape(obj.geo.coord(obj.geo.topol{obj.m+1}(obj.dom.ind(:),:),:),[obj.dom.num,obj.m+1,obj.geo.d]);
            Z=permute(Z,[1,3,2]);
            T=tensorprod(Z,n,3,1);
            T=permute(T,[1,3,2]);

            [E,G,I]=size(T);
            T=reshape(T,E*G,I);

            if obj.p<=0
                obj.coord=T;
                obj.topol=reshape(1:size(T,1),[obj.dom.num,obj.ref.num]);
            elseif obj.p>0                   
                [C,~,IC]=uniquetol(T,obj.geo.tol,'ByRows',true);
                obj.topol=reshape(IC,E,G);
                obj.coord=C;
            end
        end

        function eshaperef(obj)
            if obj.p==0
                obj.ref.coord=ones(1,obj.m)./(obj.m+1);
                obj.ref.shape={ones(obj.gauss.num,1)};
                obj.ref.num=1;
                return;
            end

            q=abs(obj.p);
            Y=1/factorial(obj.m)*prod(q+1:q+obj.m);

            if q~=0         
                args=repmat({0:q},1,obj.m);
                [grids{1:obj.m}]=ndgrid(args{:});
                ex=cell2mat(cellfun(@(x) x(:),grids,'UniformOutput',false));
                ex=ex(sum(ex,2)<=q,:);
                ex=sortrows(ex);
                obj.ref.coord=ex;
    
                V=reshape(obj.ref.coord/q,[1,Y,obj.m]).^reshape(ex,[Y,1,obj.m]);
                V=reshape(V,[],obj.m);
                V=prod(V,2);
                V=reshape(V,[Y,Y]);
                C=inv(V);
    
                alpha=reshape(ex,[Y,1,obj.m]);
                beta=reshape(ex,[1,Y,obj.m]);
                d=alpha-beta;
                d(d<0)=1;
            else
                obj.ref.coord=ones(1,obj.m)./(obj.m+1);
                obj.ref.shape{q,1}=ones(1,g);
            end

            obj.ref.num=size(obj.ref.coord,1);
            obj.ref.shape=cell(q+1,1);

            D=reshape(obj.gauss.coord,[1,1,obj.gauss.num,obj.m]).^reshape(d,[Y,Y,1,obj.m]);
            D=D.*reshape(((alpha-beta)>=0)*1,[Y,Y,1,obj.m]);
            D=D.*reshape(factorial(alpha)./factorial(d),[Y,Y,1,obj.m]);
            D=reshape(D,[],obj.m);
            D=prod(D,2);
            D=reshape(D,[Y,Y,obj.gauss.num]);
            D=tensorprod(C,D,2,1);
    
            for r=0:q
                id=(sum(obj.ref.coord,2)==r);
                refcoord0=obj.ref.coord(id,:);
                val=D(:,id,:);

                if r==0
                    obj.ref.shape{r+1}=squeeze(val)';
                else
                    tmp=zeros([Y,repmat(obj.m,1,r),obj.gauss.num]);

                    for a=1:nnz(id)
                        der=cell(obj.m,1);
                        for i=1:obj.m
                            der{i}=repmat(i,1,refcoord0(a,i));
                        end
                        der=horzcat(der{:});
                        ders=perms(der);
        
                        for der=ders'
                            i=num2cell(der');
                            tmp(:,i{:},:)=squeeze(val(:,a,:));
                        end
                    end

                    tmp=reshape(tmp,[Y,obj.m^r*obj.gauss.num]);
                    tmp=tmp';
                    tmp=reshape(tmp,[repmat(obj.m,1,r),obj.gauss.num,Y]);
                    obj.ref.shape{r+1}=tmp;
                end
            end
        end
    
        function eshape(obj)
            function out=addbasis(in)
                if obj.basis.dim==1
                    if isscalar(obj.basis.B)
                        out=obj.basis.B*in;
                    else
                        out=tensorprod(obj.basis.B,in);
                    end
                else
                    szN=size(obj.basis.B);
                    tp=reshape(in,[ones(1,ndims(obj.basis.B)-1),size(in),1]);
                    base=reshape(obj.basis.B,[szN(1:end-1),ones(1,ndims(in)),obj.basis.dim]);
                    tp=tp.*base;
                    if obj.p~=0
                        szN=size(tp);
                        tp=reshape(tp,[szN(1:end-2),prod(szN(end-1:end))]);
                    end
                    out=tp;
                end
            end

            obj.hdshape=cell(abs(obj.p)+1,1);

            A=[-ones(1,obj.m);eye(obj.m)];
            Z=reshape(obj.geo.coord(obj.geo.topol{obj.m+1}(obj.dom.ind(:),:),:),[obj.dom.num,obj.m+1,obj.geo.d]);
            Z=permute(Z,[3,2,1]);
            ZA=permute(tensorprod(Z,A,2,1),[1,3,2]);
            J=pagepinv(ZA);

            tmp=repmat(obj.ref.shape{1}(:)',[obj.dom.num,1]);
            tmp=reshape(tmp,[obj.dom.num,size(obj.ref.shape{1})]);
            obj.hdshape{1}=addbasis(tmp);

            JJ=J;

            for r=1:abs(obj.p)
                tmp=obj.ref.shape{r+1};
                tmp=tensorprod(tmp,JJ,1:r,1:r);
                tmp=permute(tmp,[3:ndims(tmp),1,2]);

                obj.hdshape{r+1}=addbasis(tmp);

                szN=size(JJ);
                szN=szN(1:end-1);
                E=obj.dom.num;
                J1=reshape(J,[ones(1,2*r),size(J,1),size(J,2),E]);
                JJ1=reshape(JJ,[szN,1,1,E]);
                JJ=JJ1.*J1;
                JJ=permute(JJ,[1:r,2*r+1,r+1:2*r,2*r+2,2*r+3]);
            end

            obj.shape=cell(obj.geo.m,abs(obj.p)+1);

            for r=1:abs(obj.p)+1
                N=obj.hdshape{r};
                szN=size(N);

                for j=1:obj.geo.m

                    if j<=obj.m
                        T=obj.geo.gauss{obj.m}.trace{j}.conn(:,1);
                    else
                        T=obj.geo.gauss{j}.trace{obj.m}.conn(:,2);
                    end

                    if obj.ref.ndof==1
                        szM=szN;
                        szM(end-1)=obj.geo.num(obj.m+1);
                        M=zeros(szM);
                        if obj.basis.dim==1
                            M(obj.dom.ind,:)=N;
                            obj.shape{j,r}=M(T)';
                        else
                            ix=repmat({':'},1,ndims(M)-2);
                            M(ix{:},obj.dom.ind,:)=N;
                            M=reshape(M,[prod(szM(1:end-2)),szM(end-1)*szM(end)]);
                            tmp=M(:,T);
                            obj.shape{j,r}=reshape(tmp,[szM(1:end-2),length(T)]);
                        end
                    else
                        szM=szN;
                        szM(end-2)=obj.geo.num(obj.m+1);
                        M=zeros(szM);
                        ix=repmat({':'},1,ndims(M)-3);
                        M(ix{:},obj.dom.ind,:,:)=N;
                        n=ndims(M);
                        M=permute(M,[1:n-3,n,n-2,n-1]);
                        M=reshape(M,[prod(szM(1:end-3))*szM(end),szM(end-2)*szM(end-1)]);
                        tmp=M(:,T);
                        obj.shape{j,r}=reshape(tmp,[szM(1:end-3),szM(end),length(T)]);
                    end
                end
            end

            nP=abs(obj.p)+1;
            om=obj.m;
            nd=obj.ref.ndof;
            for r=1:nP
                N=obj.hdshape{r};
                szN=size(N);
                if nd==1
                    szM=szN;
                    szM(end-1)=obj.geo.num(om+1);
                    M=zeros(szM);
                    if obj.basis.dim==1
                        M(obj.dom.ind,:)=N;
                        for j=1:om
                            if j<=om
                                T=obj.geo.gauss{om}.trace{j}.conn(:,1);
                            else
                                T=obj.geo.gauss{j}.trace{om}.elem(:,2);
                            end
                            obj.shape{j,r}=M(T)';
                        end
                    else
                        ix=repmat({':'},1,ndims(M)-2);
                        M(ix{:},obj.dom.ind,:)=N;
                        M=reshape(M,[prod(szM(1:end-2)),szM(end-1)*szM(end)]);
                        for j=1:om
                            T=obj.geo.gauss{om}.trace{j}.conn(:,1);
                            tmp=M(:,T);
                            obj.shape{j,r}=reshape(tmp,[szM(1:end-2),length(T)]);
                        end
                    end
                else
                    szM=szN;
                    szM(end-2)=obj.geo.num(om+1);
                    M=zeros(szM);
                    ix=repmat({':'},1,ndims(M)-3);
                    M(ix{:},obj.dom.ind,:,:)=N;
                    n=ndims(M);
                    M=permute(M,[1:n-3,n,n-2,n-1]);
                    M=reshape(M,[prod(szM(1:end-3))*szM(end),szM(end-2)*szM(end-1)]);
                    for j=1:om
                        T=obj.geo.gauss{om}.trace{j}.conn(:,1);
                        tmp=M(:,T);
                        obj.shape{j,r}=reshape(tmp,[szM(1:end-3),szM(end),length(T)]);
                    end
                end
            end
        end

        function obj=space(dom,p,options)
            arguments
                dom
                p 
                options.type="Sca"
            end

            obj.dom=dom;
            obj.geo=dom.geo;
            obj.m=dom.m;
            obj.p=p;

            obj.gauss=obj.geo.gauss{obj.m}.ref;

            if isstring(options.type)
                if strcmp(options.type,"Sca")
                    obj.basis.B=1;
                    obj.basis.dim=1;
                    obj.basis.size=1;
                elseif strcmp(options.type,"Vec")
                    obj.basis.B=eye(obj.geo.d);
                    obj.basis.dim=obj.geo.d;
                    obj.basis.size=[obj.geo.d,1];
                elseif strcmp(options.type,"Mat")
                    B=eye(obj.geo.d^2);
                    obj.basis.B=reshape(B,[obj.geo.d,obj.geo.d,obj.geo.d^2]);
                    obj.basis.dim=obj.geo.d^2;
                    obj.basis.size=[obj.geo.d,obj.geo.d];
                elseif strcmp(options.type,"Sym")
                    B=eye(obj.geo.d*obj.geo.d);
                    B=reshape(B,[obj.geo.d,obj.geo.d,obj.geo.d^2]);
                    B=(B+permute(B,[2,1,3]))/2;
                    norm=pagenorm(B,"fro");
                    B=B./norm;
                    mask=triu(true(obj.geo.d)); 
                    B=B(:,:,mask(:));
                    obj.basis.B=B;
                    obj.basis.dim=fix(obj.geo.d*(obj.geo.d+1)/2);
                    obj.basis.size=[obj.geo.d,obj.geo.d];
                elseif strcmp(options.type,"DevSym")
                    d=obj.geo.d;

                    [I,J]=find(triu(ones(d),1)); 
                    noff=length(I);

                    id1=sub2ind([d,d],I,J);
                    id2=sub2ind([d,d],J,I);
                    
                    Boff=zeros(d,d,noff);
                    Boff(id1+d*d*(0:noff-1)')=1/sqrt(2);
                    Boff(id2+d*d*(0:noff-1)')=1/sqrt(2);
                    
                    spacediag=null(ones(1,d));
                    ndiag=d-1;
                    
                    Bdiag=zeros(d,d,ndiag);
                    for k=1:ndiag
                        Bdiag(:,:,k)=diag(spacediag(:,k));
                    end

                    B=cat(3,Boff,Bdiag);

                    obj.basis.B=B;
                    obj.basis.dim=d*(d+1)/2-1;

                    obj.basis.size=[d,d];
                elseif strcmp(options.type,"Dev")
                    d=obj.geo.d;

                    Boff=reshape(eye(d^2),[d,d,d^2]); 
                    isdiag=logical(eye(d));
                    idx=~isdiag(:);
                    Boff=Boff(:,:,idx);
                    
                    [Q,~]=qr(ones(d,1)); 
                    diagvec=Q(:,2:end); 
                    
                    Bdiag=reshape(diagvec,[d,1,d-1]).*eye(d);
                    
                    B=cat(3,Boff,Bdiag);
                    
                    obj.basis.B=B;
                    obj.basis.dim=d-1;

                    obj.basis.size=[d,d];
                elseif strcmp(options.type,"Skew")
                    d=obj.geo.d;

                    [I,J]=find(triu(ones(d),1));
                    num=length(I);
                    
                    B=zeros(d,d,num);
                    
                    offsets=(d*d*(0:num-1))';
                    id=sub2ind([d,d],I,J)+offsets;
                    idt=sub2ind([d,d],J,I)+offsets;

                    val=1/sqrt(2);
                    B(id)=val;
                    B(idt)=-val;
                    
                    obj.basis.B=B;
                    obj.basis.dim=fix(d*(d-1)/2);

                    obj.basis.size=[d,d];
                end
            else
                obj.basis.size=options.type;
                obj.basis.dim=prod(options.type);
                obj.basis.B=reshape(eye(obj.basis.dim),[options.type,prod(options.type)]);
            end

            obj.eshaperef;
            obj.enode;
            obj.enum;
            obj.eind;
            obj.eshape;
        end
    end
end