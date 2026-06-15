classdef domain< handle
    properties
        geo
        m
        ind
        num
    end
    methods
        function out=sub(obj,f)
            incl=f(obj.geo.coord');
            mask=all(incl(obj.geo.topol{obj.m+1}(obj.ind,:)),2);

            out=domain(obj.geo,obj.m,obj.ind(mask));
        end

        function out=bound(obj)
            mask=sum(obj.geo.adj{obj.m}(obj.ind,:),1)==1;
            out=domain(obj.geo,obj.m-1,find(mask));
        end

        function out=skel(obj)
            mask=sum(obj.geo.adj{obj.m}(obj.ind,:),1)==2;
            out=domain(obj.geo,obj.m-1,find(mask));
        end

        function obj=domain(geo,m,ind)
            assert(m>0,"un dominio deve avere dimensione maggiore o uguale ad 1");

            obj.geo=geo;
            obj.m=m;
            obj.ind=ind;

            obj.num=length(ind);
        end
    end
end