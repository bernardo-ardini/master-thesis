classdef field < handle
    properties
        space % space
        geo % geometry
        dom % domain
        dof % dof(a) ath dof
        hdval % hdval{r}(I1,...,IN,i1,...,ir,e,g) derivative of the component (I1,...,IN) of the field along direction (i1,...,ik) computed at gauss point g of element e
        val % val{k,r}(I1,...,IN,i1,...,ir,G)
    end
    methods
        function obj=field(space) % create a field of a given space
            obj.dof=zeros(space.ndof,1);
            obj.space=space;
            obj.geo=space.geo;
            obj.dom=space.dom;
        end

        function map(obj,f) % compute the dofs of the field given f, given x(i,a) the funct handle f should return an array of size (I1,...,IN,size(x,2))
            tmp=f(obj.space.coord')';
            obj.dof=tmp(:);
        end

        function eval(obj) % update val
            obj.hdval=cell(abs(obj.space.p)+1,1);
            for r=1:abs(obj.space.p)+1
                locdof=obj.dof(obj.space.hdind);
                if obj.space.ref.ndof==1
                    szV=ones(1,ndims(obj.space.hdshape{r}));
                    szV(end-1)=obj.dom.num;
                    locdof=reshape(locdof,szV);
                    obj.hdval{r}=obj.space.hdshape{r}.*locdof;
                else
                    szV=ones(1,ndims(obj.space.hdshape{r}));
                    szV(end-2)=obj.dom.num;
                    szV(end)=obj.space.ref.ndof;
                    locdof=reshape(locdof,szV);
                    obj.hdval{r}=sum(obj.space.hdshape{r}.*locdof,length(szV));
                end
            end

            P=abs(obj.space.p)+1;
            md=obj.dom.m;
            obj.val=cell(obj.geo.m,P);
            for r=1:P
                V=obj.hdval{r};

                for j=1:obj.geo.d
                    if j<=md
                        T=obj.geo.gauss{md}.trace{j}.conn(:,1);
                    else
                        T=obj.geo.gauss{j}.trace{md}.conn(:,2);
                    end

                    if obj.space.gauss.num==1
                        V=reshape(V,[prod(obj.space.basis.size)*obj.geo.d^(r-1),obj.dom.num]);
                        W=zeros(prod(obj.space.basis.size)*obj.geo.d^(r-1),obj.geo.num(md+1));
                        W(:,obj.dom.ind)=V;
                    else
                        V=reshape(V,[prod(obj.space.basis.size)*obj.geo.d^(r-1),obj.dom.num,obj.space.gauss.num]);
                        W=zeros(prod(obj.space.basis.size)*obj.geo.d^(r-1),obj.geo.num(md+1),obj.space.gauss.num);
                        W(:,obj.dom.ind,:)=V;
                        W=reshape(W,[prod(obj.space.basis.size)*obj.geo.d^(r-1),obj.geo.num(md+1)*obj.space.gauss.num]);
                    end

                    tmp=W(:,T);
                    tmp=reshape(tmp,[obj.space.basis.size,repmat(obj.geo.d,1,r-1),length(T)]);
                    obj.val{j,r}=squeeze(tmp);
                end
            end
        end

        function plot(obj,options)
            arguments
                obj
                options.r=1
                options.prc="val"
                options.displ="none"
                options.scale=1
                options.trisurf=0
                options.LineStyle="-"
            end

            vcoord=round(obj.geo.coord/obj.geo.tol)*obj.geo.tol;
            ncoord=round(obj.space.coord/obj.geo.tol)*obj.geo.tol;
            
            [vind,nind]=ismember(vcoord,ncoord,'rows');

            va=obj.dof;
            va=reshape(va,obj.space.num,[]);
            va=va(nind(vind),:);
            va=va';
            va=reshape(va,[obj.space.basis.dim,size(obj.geo.coord,1)]);
            va=squeeze(va);

            if strcmp(options.prc,"val")
                va=va(:);
            elseif strcmp(options.prc,"norm")
                va=reshape(va,[],size(obj.geo.coord,1));
                va=squeeze(sqrt(sum(va.^2,1)));
                va=va(:);
            elseif isa(options.prc,"function_handle")
                va=options.prc(va);
                va=va(:);
            end

            coord=obj.geo.coord;

            if ~strcmp(options.displ,"none")
                C=obj.geo.gauss{2}.coord;
                va=options.displ.val{2,1};
                va=va';
                U=scatteredInterpolant(C(:,1),C(:,2),va(:,1));
                U.Values=va;
                U=F(coord(:,1),coord(:,2));
                coord=coord+options.scale*U;
            end

            if options.trisurf==0
                patch('Faces',obj.geo.topol{end},'Vertices',coord,'FaceVertexCData',va,'FaceColor','interp','LineStyle',options.LineStyle);
            else
                tr=trisurf(obj.geo.topol{end},coord(:,1),coord(:,2),va,'LineStyle',options.LineStyle);
                shading interp;
                if options.LineStyle=="-"
                    set(tr,'EdgeColor','k','LineWidth',0.5);
                end
            end

            pbaspect([1,1,1]);
            axis equal;
            colormap parula(10);
        end
    end
end