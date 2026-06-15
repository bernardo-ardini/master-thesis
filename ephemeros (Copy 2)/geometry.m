classdef geometry < handle
    properties
        d % dim spazio ambiente
        m % dim geometria
        topol % topol{k}(e,a) nodo dell'elemento e di dim k-1
        coord % coord(a,i) coord i nodo a
        num % num(k) num elementi di dim k-1
        adj % adj{k}(f,e) è 1 se l'elemento e di dim k-1 è contenuto nell'elemento f di dim k
        pres % pres{k}(e,a) è 1 se il nodo a è contenuto nell'elemento e di dim k
        conn % q un indice lineare, conn{k}(2,q) elemento di dim k-1 contenuto in conn{k}(1,q) elemento di dim k
        sign % q un indice lineare, sign{k}(q) è un segno scelto
        mass % mass{k}(e) massa (misura hausdorff k-1) dell'elemento e di dim k-1
        normal % normal{k}(i,q) coord i del vettore unitario, normale a conn{k}(2,q) elemento di dim k-1, tangente e uscente a conn{k}(1,q) elemento di dim k
        tol % tolerance for geometric computations
        gauss
        bari
        % gauss{k} gestione dei punti di gauss sugli elementi di dimensione k
        % gauss{k}.ref punti di gauss nell'elemento di riferiemento esclusi quelli che vengono da elementi di dimensione inferiore
        % gauss{k}.set.coord(g,i) coord i del punto di gauss locale g
        % gauss{k}.set.num numeri di punti di gauss locali
        % gauss{k}.set.weight(g) peso del punto di gauss locale g
        % gauss{k}.ref punti di gauss nell'elemento di riferiemento compresi quelli che vengono da elementi di dimensione inferiore
        % gauss{k}.ref.coord(g,i) coordinata del punto di gauss locale g nell'elemento di riferiemento
        % gauss{k}.ref.num numero di punti di gauss nell'elemento di riferimento
        % gauss{k}.coord(g,i) coord i del punto di gauss globale g
        % gauss{k}.elem(g) elemento a cui appartiene il punto di gauss g
        % gauss{k}.num numero di punti di gauss globali
        % gauss{k}.measure(g) prodotto della massa dell'elemento e del peso del punto di gauss
        % gauss{k}.trace{j}.conn(G,1) indice lineare dei punti di gauss nell'elemento e che sono appartengono ad elementi di dimensione j
        % gauss{k}.trace{j}.conn(G,2) indice del punto di gauss globale sugli elementi di dimensione k corrispondente all'indice lineare G
        % gauss{k}.trace{j}.num
        % gauss{k}.trace{j}.measure(G)
        % gauss{k}.trace{k-1}.q(G)
        % gauss{k}.trace{k-1}.normal(i,G)
    end
    methods
        function topology(obj,topol) % calcola la topologia
            obj.num=zeros(obj.m+1,1);
            obj.topol=cell(obj.m+1,1);
            obj.pres=cell(obj.m+1,1);
            obj.adj=cell(obj.m,1);
            obj.conn=cell(obj.m,1);
            obj.sign=cell(obj.m,1);

            for k=1:obj.m+1
                ind=nchoosek(1:size(topol,2),k);
                tmp=reshape(topol(:,ind(:)),[],size(ind,2));
                tmp=sort(tmp,2);
                tmp=unique(tmp,"rows");
                obj.topol{k}=tmp;

                obj.num(k)=size(tmp,1);

                ro=repmat((1:obj.num(k))',1,size(tmp,2));
                co=tmp;
                ro=ro(:);
                co=co(:);
                I=sparse(ro,co,1,obj.num(k),obj.num(1));
                obj.pres{k}=logical(I);
            end

            assert(all(obj.num>1,1),"ci devono essere almeno due elementi per ogni dimensione");

            for k=1:obj.m
                obj.adj{k}=(obj.pres{k+1}*obj.pres{k}'==k);
                [ind1,ind0]=find(obj.adj{k});
                obj.conn{k}=[ind1,ind0];
                tmp=-ones(size(ind0));
                [~,ia]=unique(ind0,'stable');
                tmp(ia)=1;
                obj.sign{k}=tmp;
            end
        end

        function analysis(obj,coord) % calcola massa e normali
            obj.coord=coord;
            obj.mass=cell(obj.m+1,1);
            obj.normal=cell(obj.m,1);
            obj.bari=cell(obj.m+1,1);

            for k=1:obj.m+1           
                X=reshape(obj.coord(obj.topol{k}(:),:),[obj.num(k),k,obj.d]);
                X=permute(X,[3,2,1]);
                V=X-X(:,1,:);
                V=V(:,2:end,:);
                tmp=pagemtimes(V,'transpose',V,'none');
                tmp=squeeze(prod(pageeig(tmp),1));
                obj.mass{k}=1/factorial(k-1)*sqrt(tmp);
                obj.bari{k}=reshape(mean(X,2),obj.d,[])';
            end

            for k=1:obj.m
                [ind1,ind0]=find(obj.adj{k});
                N=size(ind1,1);

                topol0=obj.topol{k}(ind0,:);
                topol1=obj.topol{k+1}(ind1,:); 

                X0=reshape(obj.coord(topol0(:),:),[N,k,obj.d]);
                X1=reshape(obj.coord(topol1(:),:),[N,k+1,obj.d]);

                X0=permute(X0,[3,2,1]);
                X1=permute(X1,[3,2,1]);

                y=sum(X1,2)-sum(X0,2);
                T=X0(:,2:end,:)-X0(:,1,:);
                v=y-X0(:,1,:);
                C=pagemtimes(T,'transpose',T,'none');
                Tv=pagemtimes(T,'transpose',v,'none');
                la=pagemldivide(C,Tv);
                tmp=pagemtimes(T,la)-v;
                tmp=tmp./pagenorm(tmp);
                tmp=squeeze(tmp);
                
                obj.normal{k}=tmp;
            end
        end

        function out=diamest(obj)
            mass0=sum(obj.mass{end},1);
            out=mass0^(1/obj.m);
        end

        function init(obj,topol,coord,options) % inizializza fornendo topol(e,a) e coord(a,i)
            arguments
                obj 
                topol 
                coord 
                options.gord=1
            end
            obj.tol=1e-7;
            obj.d=size(coord,2);
            obj.m=size(topol,2)-1;
            obj.topology(topol);
            obj.analysis(coord);
            obj.setquad(options.gord);
            obj.egauss;
        end

        function readmodel(obj,model) % legge un PDEModel di MATLAB
            obj.init(model.Mesh.Elements',model.Mesh.Nodes');
        end

        function readgmsh(obj,filename)
            system(sprintf('gmsh %s',filename));
            gmsh;
            obj.init(msh.TRIANGLES(:,1:3),msh.POS(:,1:2));
        end

        function out=whole(obj) % restituisce il dominio corrispondente all'intera geometria
            out=domain(obj,obj.m,1:obj.num(end));
        end

        function setquad(obj,gord)
            assert(obj.m==2,"funzionalità non ancora implementata per m diverso da 2");

            obj.gauss=cell(obj.m,1);

            if gord==0
                obj.gauss{1}.set.coord=0.5;
                obj.gauss{1}.set.num=1;
                obj.gauss{1}.set.weight=1;
    
                obj.gauss{2}.set.coord=[1/3,1/3];
                obj.gauss{2}.set.num=1;
                obj.gauss{2}.set.weight=1;
            elseif gord==1
                obj.gauss{1}.set.coord=[0.5-1/(2*sqrt(3));0.5+1/(2*sqrt(3))];
                obj.gauss{1}.set.num=2;
                obj.gauss{1}.set.weight=[0.5;0.5];
    
                obj.gauss{2}.set.coord=[1/3,1/3;0.6,0.2;0.2,0.6;0.2,0.2];
                obj.gauss{2}.set.num=4;
                obj.gauss{2}.set.weight=[-27/48;25/48;25/48;25/48];
            end
        end

        function egauss(obj)
            for k=1:obj.m
                Z=reshape(obj.coord(obj.topol{k+1}(:,:),:),[obj.num(k+1),k+1,obj.d]);
                Z=permute(Z,[1,3,2]);
                a=zeros(k+1,1);
                a(1)=1;
                A=[-ones(1,k);eye(k)];
                n=a+A*obj.gauss{k}.set.coord';
                T=tensorprod(Z,n,3,1);
                T=permute(T,[1,3,2]);
                [E,G,I]=size(T);
                T=reshape(T,E*G,I);
                obj.gauss{k}.coord=T;
                obj.gauss{k}.elem=repmat((1:E)',[G,1]);

                obj.gauss{k}.num=obj.num(k+1)*obj.gauss{k}.set.num;
                
                tmp=obj.mass{k+1}*obj.gauss{k}.set.weight';
                obj.gauss{k}.measure=tmp(:);

                gcoord=cell(k,1);
                gnum=zeros(k,1);
                for j=1:k
                    if j==k
                        gcoord{j}=obj.gauss{k}.set.coord;
                        gnum(j)=obj.gauss{k}.set.num;
                    else
                        gamma=cell(nchoosek(k+1,j));
                        a=zeros(j+1,1);
                        a(1)=1;
                        A=[-ones(1,j);eye(j)];
                        n=a+A*obj.gauss{j}.set.coord';
                        c=[eye(k),zeros(k,1)];
                        comb=nchoosek(1:k+1,j+1);
                        for l=1:nchoosek(k+1,j)
                            Z=c(:,comb(l,:));
                            gamma{l}=Z*n;
                        end
                        gamma=horzcat(gamma{:});
                        gamma=gamma';
                        gcoord{j}=gamma;
                        gnum(j)=size(gamma,1);
                    end
                end
                obj.gauss{k}.ref.coord=vertcat(gcoord{:});
                obj.gauss{k}.ref.num=sum(gnum,1);

                obj.gauss{k}.trace=cell(k,1);

                Z=reshape(obj.coord(obj.topol{k+1}(:,:),:),[obj.num(k+1),k+1,obj.d]);
                Z=permute(Z,[1,3,2]);

                a=zeros(k+1,1);
                a(1)=1;
                A=[-ones(1,k);eye(k)];
                n=a+A*obj.gauss{k}.ref.coord';

                T=tensorprod(Z,n,3,1);
                T=permute(T,[1,3,2]);

                [E,G,I]=size(T);
                T=reshape(T,[E*G,I]);

                obj.gauss{k}.all.coord=T;
                obj.gauss{k}.all.num=size(T,1);

                for j=1:k

                    C=obj.gauss{j}.coord;

                    % dist=sqrt(sum((T-C).^2,3));
                    % [ind0,ind1]=find(dist<obj.tol);

                    Tr=round(T/obj.tol)*obj.tol;
                    Cr=round(C/obj.tol)*obj.tol;
                    
                    % [LIA,LOCb]=ismember(Tr,Cr,'rows');
                    [LIA,LOCb]=ismember(Tr,Cr,'rows');
                    ind0=find(LIA);
                    ind1=LOCb(LIA);

                    obj.gauss{k}.trace{j}.coord=T(ind0,:);
                    % obj.gauss{k}.trace{j}.conn=[ind0,ind1];
                    obj.gauss{k}.trace{j}.conn=zeros(size(ind0,1),2);
                    obj.gauss{k}.trace{j}.conn(:,1)=ind0;
                    obj.gauss{k}.trace{j}.num=size(obj.gauss{k}.trace{j}.conn,1);
                    obj.gauss{k}.trace{j}.measure=obj.gauss{j}.measure(ind1);

                    A=repmat((1:E)',[G,1]);
                    A=A(ind0);
                    B=[A,obj.gauss{j}.elem(ind1)];
                    obj.gauss{k}.trace{j}.elem=B;

                    if j==k-1
                        [~,q]=ismember(B,obj.conn{k},'rows');
                        obj.gauss{k}.trace{j}.q=q;
                        obj.gauss{k}.trace{j}.normal=obj.normal{k}(:,obj.gauss{k}.trace{j}.q);
                        obj.gauss{k}.trace{j}.sign=obj.sign{k}(obj.gauss{k}.trace{j}.q);
                    end
                end
            end

            for j=1:obj.m
                C=obj.gauss{j}.all.coord;

                for k=j:obj.m
                    T=obj.gauss{k}.trace{j}.coord;

                    Tr=round(T/obj.tol)*obj.tol;
                    Cr=round(C/obj.tol)*obj.tol;

                    [~,ind]=ismember(Tr,Cr,'rows');

                    obj.gauss{k}.trace{j}.conn(:,2)=ind;
                end
            end
        end

        function plot(obj,options)
            arguments
                obj
                options.label=[];
                options.facecolor="g"
                options.field=[];
            end

            assert(obj.d==2);
            assert(obj.m==2);

            if isempty(options.field)
                patch('Faces',obj.topol{end},'Vertices',obj.coord,'FaceColor',options.facecolor);
            else
                patch('Faces',obj.topol{end},'Vertices',obj.coord,'FaceVertexCData',options.field.val,'FaceColor','interp');
                if options.field.max==options.field.min
                    clim([options.field.min,options.field.max+1]);
                else
                    clim([options.field.min,options.field.max]);
                end
            end
            pbaspect([1,1,1]);
            axis equal;
            hold on;

            for k=options.label
                x1=mean(reshape(obj.coord(obj.topol{k},1),size(obj.topol{k})),2);
                x2=mean(reshape(obj.coord(obj.topol{k},2),size(obj.topol{k})),2);
                text(x1,x2,string(1:obj.num(k)));
            end
        end

        function plotfield(obj,vals,k,prc,options)
            arguments
                obj
                vals
                k
                prc
                options.displ=[]
                options.scale=1
            end

            ccoord=obj.gauss{k}.trace{k}.coord;
            vals=prc(vals);
            vals=vals(:);

            x=ccoord(:,1);
            y=ccoord(:,2);

            gcoord=obj.coord;

            F=scatteredInterpolant(x,y,vals,'linear','nearest');
            vals=F(gcoord(:,1),gcoord(:,2));

            if ~isempty(options.displ)
                uvcoord=round(obj.coord/obj.tol)*obj.tol;
                uncoord=round(options.displ.coord/obj.tol)*obj.tol;

                [uvind,unind]=ismember(uvcoord,uncoord,'rows');

                uva=options.displ.dof;
                uva=reshape(uva,options.displ.num,[]);
                uva=uva(unind(uvind),:);
                uva=uva';
                uva=reshape(uva,[options.displ.dim,size(obj.coord,1)]);
                uva=squeeze(uva);

                gcoord=gcoord+(options.scale*uva)';
            end

            trisurf(obj.topol{end},gcoord(:,1),gcoord(:,2),vals);
            view(2);
        end
    end
end