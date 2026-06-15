classdef geometryold < handle
    properties
        d
        m
        coord
        topol
        connect
        num
        
        sets

        normal
        mass
    end
    methods
        function tidy(obj)
            obj.num(obj.m+1)=size(obj.topol{obj.m+1},1);

            for k=1:obj.m+1
                tmp=arrayfun(@(e) nchoosek(obj.topol{obj.m+1}(e,:),k),1:obj.num(obj.m+1),'UniformOutput',false);
                tmp=vertcat(tmp{:});
                tmp=sort(tmp,2);
                tmp=unique(tmp,"rows");
                obj.topol{k}=tmp;
                obj.num(k)=size(obj.topol{k},1);
            end

            for k=1:obj.m
                tmp=arrayfun(@(b) sparse(all(ismember(obj.topol{k},obj.topol{k+1}(b,:)),2))',1:obj.num(k+1),'UniformOutput',false);
                obj.connect{k}=vertcat(tmp{:});
            end
        end

        function out=diamest(obj)
            mass0=sum(obj.mass{end},1);
            out=mass0^(1/obj.m);
        end

        function compute(obj)
            if obj.m==0
                return;
            end

            obj.mass=cell(obj.m,1);
            obj.normal=cell(obj.m,1);

            for k=1:obj.m
                obj.mass{k}=zeros(obj.num(k+1),1);

                for a=1:obj.num(k+1)
                    X=obj.coord(obj.topol{k+1}(a,:),:)';
                    V=X-X(:,1);
                    V=V(:,2:end);
                    obj.mass{k}(a)=1/factorial(k)*sqrt(det(V'*V));
                end

                obj.normal{k}=cell(obj.num(k),1);

                for a=1:obj.num(k)
                    ind=1:obj.num(k+1);
                    ind=ind(obj.connect{k}(:,a));

                    i=0;
                    obj.normal{k}{a}=cell(length(ind),1);

                    for b=ind
                        i=i+1;

                        n=obj.topol{k+1}(b,~ismember(obj.topol{k+1}(b,:),obj.topol{k}(a,:)));

                        x=obj.coord(n,:)';
                        y=obj.coord(obj.topol{k}(a,1),:)';
                        v=y-x;

                        if k>1
                            Y=obj.coord(obj.topol{k}(a,2:end),:)'-y;
                            N=(eye(obj.d)-Y*((Y'*Y)\Y'));
                            v=N*v;
                        end
                        
                        obj.normal{k}{a}{i}=v/norm(v);
                    end

                    if isscalar(obj.normal{k}{a})
                        obj.normal{k}{a}=obj.normal{k}{a}{1};
                    end
                end
            end
        end

        function init(obj,topol,coord)
            obj.d=size(coord,2);
            obj.m=size(topol,2)-1;
            obj.topol{obj.m+1}=topol;
            obj.coord=coord;
            obj.tidy;
            obj.compute;
            set.m=obj.m;
            set.ind{obj.m+1}=1:obj.num(end);
            obj.addset(set);
        end

        function out=addset(obj,set)
            set.ind{end}=reshape(set.ind{end},1,[]);
            set.num(set.m+1)=length(set.ind{end});
            set.topol{set.m+1}=obj.topol{set.m+1}(set.ind{set.m+1},:);
            node=1:obj.num(1);
            mask=ismember(node,obj.topol{set.m+1}(set.ind{end},:));
            node=node(mask);
            for k=1:set.m
                set.ind{k}=1:obj.num(k);
                mask=all(ismember(obj.topol{k},node),2);
                set.ind{k}=set.ind{k}(mask);
                set.num(k)=length(set.ind{k});
                set.topol{k}=obj.topol{k}(set.ind{k},:);
            end
            for k=1:set.m
                set.connect{k}=obj.connect{k}(set.ind{k+1},set.ind{k});
                set.mass{k}=obj.mass{k}(set.ind{k+1});
            end
            obj.sets{end+1}=set;
            out=length(obj.sets);
        end

        function readgmsh(obj,filename)
            system(sprintf('gmsh %s',filename));
            gmsh;
            obj.init(msh.TRIANGLES(:,1:3),msh.POS(:,1:2));
        end

        function readmodel(obj,model)
            obj.init(model.Mesh.Elements',model.Mesh.Nodes');
        end

        function out=isin(obj,x)
            assert(obj.m==obj.d);
            assert(obj.d==2);

            [~,i]=min(vecnorm(-obj.coord'+x));
            el=any(obj.topol{end}==i,2);

            for e=find(el)'
                Z=obj.coord(obj.topol{end}(e,:),:)';
                if inpolygon(x(1),x(2),Z(1,:),Z(2,:))
                    out=true;
                    return;
                end
            end

            out=false;
        end

        function out=extract(obj,q,f)
            % find nodes
            b=zeros(obj.num(1),1);
            for a=1:obj.num(1)
                b(a)=f(obj.coord(a,:));
            end
            nod=1:1:obj.num(1);
            nod=nod(b==1);

            % find elements of dom that has node in nod
            set=obj.sets{q};
            mask=all(ismember(set.topol{end},nod),2);

            % create subdomain
            sub.m=set.m;
            sub.ind{sub.m+1}=set.ind{end}(mask);

            % add subdomain
            out=obj.addset(sub);
        end

        function out=bound(obj,q)
            set=obj.sets{q};
            mask=sum(set.connect{end},1)==1;

            sub.m=set.m-1;
            sub.ind{sub.m+1}=set.ind{end-1}(mask);

            out=obj.addset(sub);
        end

        function out=int(obj,q)
            set=obj.sets{q};
            mask=sum(set.connect{end},1)>1;

            sub.m=set.m-1;
            sub.ind{sub.m+1}=set.ind{end-1}(mask);

            out=obj.addset(sub);
        end

        function plot(obj,options)
            arguments
                obj
                options.node=false;
                options.element=false;
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

            if options.node
                text(obj.coord(:,1),obj.coord(:,2),string(1:obj.num(1)));
            end

            if options.element
                x1=mean(reshape(obj.coord(obj.topol{end},1),size(obj.topol{end})),2);
                x2=mean(reshape(obj.coord(obj.topol{end},2),size(obj.topol{end})),2);
                text(x1,x2,string(1:obj.num(end)));
            end
        end
    end
end