classdef multi <handle
    properties
        num
        spc
        scale
        rescomp
        stiffblk
    end
    methods
        function add(obj,spc,options)
            arguments
                obj 
                spc 
                options.scale=[1,1]
            end
            obj.num=obj.num+1;
            obj.spc{end+1}=spc;
            obj.rescomp{end+1}=[];
            obj.stiffblk{end+1,end+1}=[];
            obj.scale(end+1,:)=options.scale;
        end

        function obj=multi()
            obj.num=0;
            obj.spc={};
            obj.scale=[];
            obj.rescomp={};
            obj.stiffblk={};
        end

        function out=res(obj)
            tmp=cell(obj.num,1);
            for i=1:length(obj.spc)
                if isempty(obj.rescomp{i})
                    tmp{i}=zeros(obj.spc{i}.ndof,1);
                else
                    tmp{i}=1/obj.scale(i,1)*obj.rescomp{i};
                end
            end
            out=vertcat(tmp{:});
        end

        function out=stiff(obj)
            tmp=cell(obj.num,obj.num);
            for i=1:length(obj.spc)
                for j=1:length(obj.spc)
                    if isempty(obj.stiffblk{i,j})
                        tmp{i,j}=sparse(obj.spc{i}.ndof,obj.spc{j}.ndof);
                    else
                        tmp{i,j}=obj.scale(j,2)/obj.scale(i,1)*obj.stiffblk{i,j};
                    end
                end
            end
            out=cell2mat(tmp);
        end

        function out=comp(obj,vec,i)
            dims=arrayfun(@(j) obj.spc{j}.ndof,1:obj.num);
            start=sum(dims(1:i-1))+1;
            stop=sum(dims(1:i));
            ind=start:stop;
            out=obj.scale(i,2)*vec(ind);
        end

        function out=blk(obj,mat,i,j)
            dims=arrayfun(@(j) obj.spc{j}.ndof,1:obj.num);
            start=sum(dims(1:i-1))+1;
            stop=sum(dims(1:i));
            row=start:stop;
            start=sum(dims(1:j-1))+1;
            stop=sum(dims(1:j));
            col=start:stop;
            out=obj.scale(i,2)/obj.scale(i,1)*mat(row,col);
        end
    end
end