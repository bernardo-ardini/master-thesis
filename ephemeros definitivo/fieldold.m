classdef fieldold < handle
    properties
        geo
        doms
        fe

        dom0
        set0

        gauss
        node

        dof
        numdofloc
        numdof
        ind
        elem

        shapeloc
        shape

        hist
        proc        
    end

    methods
        function init(obj)
            obj.dom0=obj.doms{1};
            obj.set0=obj.geo.sets{obj.dom0.q};
            obj.computeshapeloc;
            obj.findgauss;

            m=obj.set0.m;

            if strcmp(obj.fe.space,"DevSym")
                obj.fe.size=fix((obj.geo.d+1)*obj.geo.d/2-1);
            else
                obj.fe.size=obj.fe.space;
            end

            if obj.fe.p==0
                obj.node.N=obj.set0.num(end);
                obj.numdofloc=size(obj.shapeloc{1,1},1);
            elseif obj.fe.p>=1
                p=abs(obj.fe.p);
                obj.numdofloc=1/factorial(m)*prod(p+1:p+m);
                obj.node.N=0;
                for k=0:m
                    if k>obj.fe.p-1
                        I=0;
                    else
                        I=nchoosek(obj.fe.p-1,k);
                    end
                    obj.node.num(k+1)=I;
                    obj.node.N=obj.node.N+obj.set0.num(k+1)*I;
                end
            elseif obj.fe.p<=-1
                p=abs(obj.fe.p);
                obj.numdofloc=1/factorial(m)*prod(p+1:p+m);
                obj.node.N=obj.numdofloc*obj.set0.num(end);
            end

            obj.numdofloc=prod(obj.fe.size)*obj.numdofloc;
            obj.numdof=prod(obj.fe.size)*obj.node.N;
            obj.dof=zeros(obj.numdof,1);

            obj.computeshape;

            obj.findnode;
            obj.computeind;

            obj.inithist;
        end

        function computeshapeloc(obj)
            p=abs(obj.fe.p);
            m=obj.geo.sets{obj.doms{1}.q}.m;
            Y=1/factorial(m)*prod(p+1:p+m);

            if p~=0             
                args=repmat({0:p},1,m);
                [grids{1:m}]=ndgrid(args{:});
                ex=cell2mat(cellfun(@(x) x(:),grids,'UniformOutput',false));
                ex=ex(sum(ex,2)<=p,:);
                ex=sortrows(ex);
                obj.node.ex=ex;
    
                V=reshape(ex/p,[1,Y,m]).^reshape(ex,[Y,1,m]);
                V=reshape(V,[],m);
                V=prod(V,2);
                V=reshape(V,[Y,Y]);
                C=inv(V);
    
                alpha=reshape(ex,[Y,1,m]);
                beta=reshape(ex,[1,Y,m]);
                d=alpha-beta;
                d(d<0)=1;
            else
                obj.node.ex=ones(1,m)./(m+1);
            end

            obj.shapeloc=cell(length(obj.doms),p+1);

            for q=1:length(obj.doms)
                dom=obj.geo.sets{obj.doms{q}.q};
                gau=obj.doms{q}.gauss;
                
                if dom.m==m
                    gamma=gau.coord;
                    g=gau.num;
                else
                    gamma=cell(nchoosek(m+1,dom.m));
                    a=zeros(dom.m+1,1);
                    a(1)=1;
                    A=[-ones(1,dom.m);eye(dom.m)];
                    n=a+A*gau.coord';
                    c=[eye(m),zeros(m,1)];
                    comb=nchoosek(1:m+1,dom.m+1);
                    for l=1:nchoosek(m+1,dom.m)
                        Z=c(:,comb(l,:));
                        gamma{l}=Z*n;
                    end
                    gamma=horzcat(gamma{:});
                    gamma=gamma';
                    g=size(gamma,1);
                end

                obj.gauss{q}.gamma=gamma;

                if p==0
                    obj.shapeloc{q,1}=ones(1,g);
                    continue;
                end

                D=reshape(gamma,[1,1,g,m]).^reshape(d,[Y,Y,1,m]);
                D=D.*reshape(((alpha-beta)>=0)*1,[Y,Y,1,m]);
                D=D.*reshape(factorial(alpha)./factorial(d),[Y,Y,1,m]);
                D=reshape(D,[],m);
                D=prod(D,2);
                D=reshape(D,[Y,Y,g]);
                D=tensorprod(C,D,2,1);
    
                for k=0:p
                    id=(sum(ex,2)==k);
                    ex0=ex(id,:);
                    val=D(:,id,:);
    
                    if k==0
                        obj.shapeloc{q,k+1}=val;
                    else
                        obj.shapeloc{q,k+1}=zeros([Y,repmat(m,1,k),g]);
    
                        for a=1:nnz(id)
                            der=cell(m,1);
                            for i=1:m
                                der{i}=repmat(i,1,ex0(a,i));
                            end
                            der=horzcat(der{:});
                            ders=perms(der);
            
                            for der=ders'
                                i=num2cell(der');
                                obj.shapeloc{q,k+1}(:,i{:},:)=squeeze(val(:,a,:));
                            end
                        end
                    end
                end
            end
        end

        function computeshape(obj)
            p=abs(obj.fe.p);

            m=obj.set0.m;

            obj.shape=cell(length(obj.doms),p+1);

            for q=1:length(obj.doms)
                for k=0:p
                end
            end

            for q=1:length(obj.doms)
                dom=obj.doms{q};
                for e=1:obj.geo.sets{dom.q}.num(end)
                    for g=1:dom.gauss.num
                        for k=0:p
                            fnd=squeeze(find(any(obj.gauss{q}.connect(e,g,:,:),4)));
                            for f=1:length(fnd)
                                A=[-ones(1,m);eye(m)];
                                Z=obj.geo.coord(obj.set0.topol{end}(fnd(f),:),:)';
                                ZA=Z*A;
                                J=((ZA'*ZA)\ZA');

                                M0=obj.shapeloc{q,k+1};
                                sz=size(M0);
                                ix=repmat({':'},1,1+k);
                                M0=M0(ix{:},squeeze(obj.gauss{q}.connect(e,g,fnd(f),:)));
        
                                if k>0
                                    for l=1:k
                                        M0=reshape(M0,[sz(1),repmat(m,1,k-l+1),repmat(obj.geo.d,1,l-1)]);
                                        d=ndims(M0)-l+1;
                                        M0=tensorprod(M0,J,d,1);
                                    end
                                end
        
                                M0=tensorprod(eye(prod(obj.fe.size)),M0);
                                M0=reshape(M0,[obj.fe.size,obj.numdofloc,repmat(obj.geo.d,1,k)]);
                                if k>0
                                    n=ndims(M0);
                                    order=[1:n-k-1,n-k+1:n,n-k];
                                    M0=permute(M0,order);
                                end
        
                                if strcmp(obj.fe.space,"DevSym")
                                    d=obj.geo.d;
                                    m=fix(d*(d+1)/2-1);
                                    L=zeros(d,d,m);
                                    
                                    if d==2
                                        L(:,:,1)=[1,0;0,-1];
                                        L(:,:,2)=[0,1;1,0];
                                    elseif d==3
                                        L(:,:,1)=[1,0,0;0,0,0;0,0,-1];
                                        L(:,:,2)=[0,0,0;0,1,0;0,0,-1];
                                        L(:,:,3)=[0,1,0;1,0,0;0,0,0];
                                        L(:,:,4)=[0,0,1;0,0,0;1,0,0];
                                        L(:,:,5)=[0,0,0;0,0,1;0,1,0];
                                    end

                                    M0=tensorprod(L,M0,3,1);
                                end

                                obj.shape{q,k+1}{f,e,g}=M0;
                            end
                        end
                    end
                end
            end
        end

        % indices

        function computeind(obj)           
            obj.ind=zeros(obj.numdofloc,obj.set0.num(end));
            for e=1:obj.set0.num(end)
                if obj.fe.p==0
                    tmp=1:obj.numdof;
                    tmp=reshape(tmp,obj.numdofloc,[]);
                    obj.ind(:,e)=tmp(:,e);
                elseif obj.fe.p>=1
                    tmp=1:obj.numdof;
                    tmp=reshape(tmp,[],obj.node.N);
                    nd=obj.node.topol(e,:);
                    tmp=tmp(:,nd);
                    tmp=tmp(:);
                    obj.ind(:,e)=tmp;
                elseif obj.fe.p<=-1
                    tmp=1:obj.numdof;
                    tmp=reshape(tmp,[obj.numdofloc,obj.set0.num(end)]);
                    tmp=tmp(:,e);
                    obj.ind(:,e)=tmp(:);
                end
            end

            obj.elem=cell(length(obj.doms),1);
            for q=1:length(obj.doms)
                set=obj.geo.sets{obj.doms{q}.q};
                obj.elem{q}=cell(set.num(end),obj.doms{q}.gauss.num);
                for e=1:set.num(end)
                    for g=1:obj.doms{q}.gauss.num
                        obj.elem{q}(e,g)=find(squeeze(any(obj.gauss{q}.connect(e,g,:,:),4)));
                    end
                end
            end
        end

        % dof

        function out=localdof(obj,q,e,g,ff)
            f=obj.elem{q}{e,g};
            f=f(ff);
            out=obj.dof;
            out=out(obj.ind(:,f));
        end

        function map(obj,f)
            if obj.fe.p==0
                out=obj.dof;
                out=reshape(out,obj.numdofloc,[]);
                for e=1:obj.set0.num(end)
                    Z=obj.geo.coord(obj.set0.topol{end}(e,:),:);
                    x=sum(Z,1)'/(obj.set0.m+1);
                    u=f(x);
                    u=u(:);
                    out(:,e)=u;
                end
                obj.dof=out(:);
            elseif obj.fe.p>=1
                out=obj.dof;
                out=reshape(out,[],obj.node.N);
                for a=1:obj.node.N
                    x=obj.node.coord(a,:)';
                    u=f(x);
                    out(:,a)=u(:);
                end
                obj.dof=out(:);
            elseif obj.fe.p<=-1
                out=zeros(prod(obj.fe.size),size(obj.node.ex,1),obj.set0.num(end));

                for ee=1:obj.set0.num(end)
                    e=obj.set0.ind{end}(ee);
                    a=zeros(obj.set0.m+1,1);
                    a(1)=1;
                    A=[-ones(1,obj.set0.m);eye(obj.set0.m)];
                    n=a+A*obj.node.ex'/abs(obj.fe.p);
                    Z=obj.geo.coord(obj.geo.topol{obj.set0.m+1}(e,:),:)';
                    Z=Z*n;
                    for a=1:size(obj.node.ex,1)
                        x=Z(:,a);
                        u=f(x);
                        u=u(:);
                        out(:,a,ee)=u;
                    end
                end

                out=out(:);
                obj.dof=out;
            end

        end

        % evalutaion

        function findgauss(obj)
            for q=1:length(obj.doms)
                set=obj.geo.sets{obj.doms{q}.q};
                gau=obj.doms{q}.gauss;

                obj.gauss{q}.coord=zeros(set.num(end),gau.num,obj.geo.d);

                mat1=zeros(set.num(end),gau.num,obj.geo.d);

                a=zeros(set.m+1,1);
                a(1)=1;
                A=[-ones(1,set.m);eye(set.m)];

                for e=1:set.num(end)
                    Z=obj.geo.coord(obj.geo.topol{set.m+1}(set.ind{end}(e),:),:)';
                    for g=1:gau.num
                        n=a+A*gau.coord(g,:)';
                        mat1(e,g,:)=Z*n;
                    end
                end

                gamma=obj.gauss{q}.gamma;
                mat2=zeros(obj.set0.num(end),size(gamma,1),obj.geo.d);

                a=zeros(obj.set0.m+1,1);
                a(1)=1;
                A=[-ones(1,obj.set0.m);eye(obj.set0.m)];

                for e=1:obj.set0.num(end)
                    n=a+A*gamma';
                    Z=obj.geo.coord(obj.set0.topol{end}(e,:),:)';
                    mat2(e,:,:)=(Z*n)';
                end

                obj.gauss{q}.coord=mat1;

                sz1=size(mat1);
                sz2=size(mat2);
                mat1=reshape(mat1,[sz1(1:2),1,1,sz1(3)]);
                mat2=reshape(mat2,[1,1,sz2]);
                obj.gauss{q}.connect=all(abs(mat1-mat2)<1e-7,5);
            end 
        end
        
        function out=eval(obj,q,k,e,g)
            m=obj.shape{q,k};
            out=cell(size(m,1),1);

            for ff=1:size(m,1)
                u=obj.localdof(q,e,g,ff);
                m=obj.shape{q,k}{ff,e,g};
                out{ff}=squeeze(tensorprod(m,u,ndims(m),1));
            end

            if isscalar(out)
                out=out{1};
            end
        end

        function out=jump(obj,q,k,e,g)
            ev=obj.eval(q,k,e,g);
            assert(iscell(ev) && length(ev)==2);
            out=ev{2}-ev{1};
        end

        function out=mean(obj,q,k,e,g)
            ev=obj.eval(q,k,e,g);
            assert(iscell(ev));
            out=sum(cat(3,ev{:}),3)/length(ev);
        end

        function findnode(obj)
            a=zeros(obj.set0.m+1,1);
            a(1)=1;
            A=[-ones(1,obj.set0.m);eye(obj.set0.m)];
            n=a+A*obj.node.ex'/abs(obj.fe.p);

            T=zeros(obj.set0.num(end),size(obj.node.ex,1),obj.geo.d);

            for e=1:obj.set0.num(end)
                Z=obj.geo.coord(obj.set0.topol{end}(e,:),:)';
                T(e,:,:)=(Z*n)';
            end

            [E,G,I]=size(T);
            T=reshape(T,E*G,I);
                
            [C,~,IC]=uniquetol(T,1e-7,'ByRows',true);
            obj.node.topol=reshape(IC,E,G);
            obj.node.coord=C;
        end

        % plot

        function out=interp(obj)
            if isempty(obj.proc)
                out=[];
                return
            end

            f0=zeros(obj.set0.num(end),obj.dom0.gauss.num,1);
            for e=1:obj.set0.num(end)
                for g=1:obj.dom0.gauss.num
                    u=obj.eval(1,1,e,g);
                    du=obj.eval(1,2,e,g);
                    f0(e,g)=obj.proc.f(u,du);
                end
            end
            f0=f0(:);
            x=reshape(obj.gauss{1}.coord,[],2);

            out=scatteredInterpolant(x,f0,'linear');
            out=out(obj.geo.coord);
        end

        % hist

        function [umin,umax,valmin,valmax]=minmax(obj)
            u=reshape(obj.dof,prod(obj.fe.size),[]);
            umin=min(u,[],2);
            umax=max(u,[],2);
            val=obj.interp;
            valmin=min(val);
            valmax=max(val);
            
            if ~isscalar(obj.fe.size)
                umin=reshape(umin,obj.fe.size);
                umax=reshape(umax,obj.fe.size);
            end
        end

        function addhist(obj,t)
            obj.hist.proc.val{end+1}=obj.interp;
            obj.hist.dof{end+1}=obj.dof;
            obj.hist.t(end+1)=t;
            [umin,umax,valmin,valmax]=obj.minmax;
            obj.hist.min=min([obj.hist.min(:),umin(:)],[],2);
            obj.hist.max=min([obj.hist.max(:),umax(:)],[],2);
            if ~isscalar(obj.fe.size)
                obj.hist.min=reshape(obj.hist.min,obj.fe.size);
                obj.hist.max=reshape(obj.hist.max,obj.fe.size);
            end
            obj.hist.proc.min=min([obj.hist.proc.min,valmin],[],2);
            obj.hist.proc.max=max([obj.hist.proc.max,valmax],[],2);
        end

        function inithist(obj)
            obj.hist.proc.val={};
            obj.hist.dof={};
            obj.hist.t=[];
            [umin,umax,valmin,valmax]=obj.minmax;
            obj.hist.min=umin;
            obj.hist.max=umax;
            obj.hist.proc.min=valmin;
            obj.hist.proc.max=valmax;
        end

        function loadhist(obj,t)
            if t<Inf
                [~,n]=min(abs(t-obj.hist.t));
                obj.dof=obj.hist.dof{n};
                obj.proc.val=obj.hist.proc.val{n};
                obj.proc.min=obj.hist.proc.min;
                obj.proc.max=obj.hist.proc.max;
            end
        end
    end
end