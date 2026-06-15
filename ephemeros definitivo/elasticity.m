classdef elasticity < handle
    properties
        geo
        GammaD
        GammaC
        OmegaM
        GammaMC

        unit

        chiM
        chi
        chibar
        v
        r
        kappa

        hist

        gamma
        dgamma

        P
        dP

        path
        setting

        f

        t
        T
    end
    methods
        % initialize, save in files and manage external actions

        function init(obj)
            dOmega.q=1;
            dOmega.gauss.coord=[1/3,1/3];
            dOmega.gauss.weights=1;
            dOmega.gauss.num=1;

            dOmega.q=1;
            dOmega.gauss.coord=[1/3,1/3];
            dOmega.gauss.weights=1;
            dOmega.gauss.num=1;

            obj.chi=field();
            obj.chi.geo=obj.geo;
            obj.chi.fe.p=1;
            obj.chi.fe.size=2;
            obj.chi.dom{1}=gauss;
            obj.chi.proc.f=obj.hist.f;
            obj.chi.init;
            obj.chi.dof=obj.geo.coord';
            obj.chi.dof=obj.chi.dof(:);
            obj.chi.inithist;

            obj.v=field();
            obj.v.Omega=obj.geo;
            obj.v.Gamma=obj.GammaD;
            obj.v.fe=fe;
            obj.v.gauss=gauss;
            obj.v.init;
            obj.v.inithist;

            gauss.coord=[1/2-sqrt(1/3)/2;1/2+sqrt(1/3)/2];
            gauss.num=2;
            gauss.weight=[0.5;0.5];

            fe.p=0;
            fe.size=2;
            obj.r=field();
            obj.r.Omega=obj.GammaD;
            obj.r.fe=fe;
            obj.r.gauss=gauss;
            obj.r.init;
            obj.r.inithist;

            fe.p=0;
            fe.size=1;
            obj.kappa=field();
            obj.kappa.Omega=obj.GammaC;
            obj.kappa.fe=fe;
            obj.kappa.gauss=gauss;
            obj.kappa.init;
            obj.kappa.inithist;

            obj.t=0;
            obj.setting.step.dt=Inf;
            obj.setting.step.lev=0;
            obj.setting.step.step=0;

            obj.unit.length=obj.geo.diamest;
            obj.unit.pressure=mean(diag(reshape(obj.dP(eye(2)),4,4)),"all");
            obj.unit.force=obj.unit.length*obj.unit.pressure;
            obj.unit.time=obj.T;
            obj.unit.velocity=obj.unit.length/obj.unit.time;

            obj.setting.FB.cx=1/obj.unit.length;
            obj.setting.FB.cy=1/obj.unit.pressure;
        end

        function write(obj)
            obj.v.addhist(obj.t);
            obj.chi.addhist(obj.t);
            obj.r.addhist(obj.t);
            obj.kappa.addhist(obj.t);
        end

        % residual, stiffness and solver for viscoplasticity 

        function out=residual(obj)
            row1a=zeros(obj.geo.num(end),obj.chi.gauss.num,obj.v.ndof);
            val1a=zeros(obj.geo.num(end),obj.chi.gauss.num,obj.v.ndof);
            for e=1:obj.geo.num(end)
                for g=1:obj.chi.gauss.num
                    X=obj.chi.x(e,g);
                    F=obj.chi.deval(e,g);
                    P0=obj.P(F);

                    int=squeeze(tensorprod(P0,obj.v.dM{e,g},[1,2],[1,2]));
                    int=int(:);

                    ext=obj.f(X,obj.t)'*obj.v.M{e,g};
                    ext=ext(:);

                    row1a(e,g,:)=obj.v.ind(e);
                    val1a(e,g,:)=obj.chi.gauss.weight(g)*(int-ext)*obj.geo.mass{end}(e);
                end
            end
            out1a=accumarray(row1a(:),val1a(:));

            bchi=obj.chi.trace(obj.GammaD,obj.r.gauss);
            row1b=zeros(obj.GammaD.num(end),obj.r.gauss.num,bchi.ndof);
            val1b=zeros(obj.GammaD.num(end),obj.r.gauss.num,bchi.ndof);
            row2=zeros(obj.GammaD.num(end),obj.r.gauss.num,obj.r.ndof);
            val2=zeros(obj.GammaD.num(end),obj.r.gauss.num,obj.r.ndof);
            for e=1:obj.GammaD.num(end)
                for g=1:obj.r.gauss.num
                    X=obj.r.x(e,g);
                    x0=bchi.eval(e,g);
                    xbar0=obj.chibar(X,obj.t);
                    r0=obj.r.eval(e,g);

                    row1b(e,g,:)=bchi.ind(e);
                    val1b(e,g,:)=-obj.r.gauss.weight(g)*r0'*bchi.M{e,g}*obj.GammaD.mass{end}(e);

                    row2(e,g,:)=obj.r.ind(e);
                    val2(e,g,:)=obj.r.gauss.weight(g)*(x0-xbar0)'*obj.r.M{e,g}*obj.GammaD.mass{end}(e);
                end
            end
            out1b=accumarray(row1b(:),val1b(:));
            out1b=bchi.tr'*out1b;
            out2=accumarray(row2(:),val2(:));

            bchi=obj.chi.trace(obj.GammaC,obj.kappa.gauss);
            row1c=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,bchi.ndof);
            val1c=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,bchi.ndof);
            row3=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof);
            val3=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof);
            for e=1:obj.GammaC.num(end)
                for g=1:obj.kappa.gauss.num
                    gamma0=obj.gamma(e,g);
                    dgamma0=squeeze(obj.dgamma(e,g,:));
                    dgamma0=dgamma0(:);
                    kappa0=obj.kappa.eval(e,g);

                    row1c(e,g,:)=bchi.ind(e);
                    val1c(e,g,:)=-obj.kappa.gauss.weight(g)*kappa0*dgamma0'*bchi.M{e,g}*obj.GammaC.mass{end}(e);

                    row3(e,g,:)=obj.kappa.ind(e);
                    val3(e,g,:)=obj.kappa.gauss.weight(g)*gamma0*obj.kappa.M{e,g}'*obj.GammaC.mass{end}(e);
                end
            end
            out1c=accumarray(row1c(:),val1c(:));
            out1c=bchi.tr'*out1c;
            gammadof=accumarray(row3(:),val3(:));

            out3=obj.FB(gammadof,obj.kappa.dof);

            out1=out1a+out1b+out1c;

            out1=out1/obj.unit.force;
            out2=out2/obj.unit.force;

            out=[out1;out2;out3];
        end

        function out=stiffness(obj)
            row11a=zeros(obj.geo.num(end),obj.chi.gauss.num,obj.v.ndof,obj.v.ndof);
            col11a=zeros(obj.geo.num(end),obj.chi.gauss.num,obj.v.ndof,obj.v.ndof);
            val11a=zeros(obj.geo.num(end),obj.chi.gauss.num,obj.v.ndof,obj.v.ndof);
            for e=1:obj.geo.num(end)
                for g=1:obj.chi.gauss.num
                    F=obj.chi.deval(e,g);
                    dP0=obj.dP(F);

                    int=squeeze(tensorprod(dP0,obj.v.dM{e,g},[3,4],[1,2]));
                    int=tensorprod(obj.v.dM{e,g},int,[1,2],[1,2]);

                    [ro,co]=meshgrid(obj.v.ind(e));

                    row11a(e,g,:,:)=ro;
                    col11a(e,g,:,:)=co;
                    val11a(e,g,:,:)=obj.chi.gauss.weight(g)*obj.dt*int*obj.geo.mass{end}(e);
                end
            end
            out11=sparse(row11a(:),col11a(:),val11a(:));

            bchi=obj.chi.trace(obj.GammaD,obj.r.gauss);
            row12=zeros(obj.GammaD.num(end),obj.r.gauss.num,bchi.ndof,obj.r.ndof);
            col12=zeros(obj.GammaD.num(end),obj.r.gauss.num,bchi.ndof,obj.r.ndof);
            val12=zeros(obj.GammaD.num(end),obj.r.gauss.num,bchi.ndof,obj.r.ndof);
            row21=zeros(obj.GammaD.num(end),obj.r.gauss.num,obj.r.ndof,bchi.ndof);
            col21=zeros(obj.GammaD.num(end),obj.r.gauss.num,obj.r.ndof,bchi.ndof);
            val21=zeros(obj.GammaD.num(end),obj.r.gauss.num,obj.r.ndof,bchi.ndof);
            for e=1:obj.GammaD.num(end)
                for g=1:obj.r.gauss.num
                    [ro,co]=meshgrid(bchi.ind(e),obj.r.ind(e));
                    row12(e,g,:,:)=ro';
                    col12(e,g,:,:)=co';
                    val12(e,g,:,:)=-obj.r.gauss.weight(g)*bchi.M{e,g}'*obj.r.M{e,g}*obj.GammaD.mass{end}(e);

                    [ro,co]=meshgrid(obj.r.ind(e),bchi.ind(e));
                    row21(e,g,:,:)=ro';
                    col21(e,g,:,:)=co';
                    val21(e,g,:,:)=obj.r.gauss.weight(g)*obj.dt*obj.r.M{e,g}'*bchi.M{e,g}*obj.GammaD.mass{end}(e);
                end
            end
            out12=sparse(row12(:),col12(:),val12(:));
            out12=bchi.tr'*out12;
            out21=sparse(row21(:),col21(:),val21(:));
            out21=out21*bchi.tr;

            bchi=obj.chi.trace(obj.GammaC,obj.kappa.gauss);
            row13=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,bchi.ndof,obj.kappa.ndof);
            col13=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,bchi.ndof,obj.kappa.ndof);
            val13=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,bchi.ndof,obj.kappa.ndof);
            row3=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof);
            val3=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof);
            row31=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof,bchi.ndof);
            col31=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof,bchi.ndof);
            val31=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,obj.kappa.ndof,bchi.ndof);
            for e=1:obj.GammaC.num(end)
                for g=1:obj.kappa.gauss.num
                    gamma0=obj.gamma(e,g);
                    dgamma0=squeeze(obj.dgamma(e,g,:));
                    dgamma0=dgamma0(:);

                    row3(e,g,:)=obj.kappa.ind(e);
                    val3(e,g,:)=obj.kappa.gauss.weight(g)*gamma0*obj.kappa.M{e,g}'*obj.GammaC.mass{end}(e);

                    [ro,co]=meshgrid(bchi.ind(e),obj.kappa.ind(e));
                    row13(e,g,:,:)=ro';
                    col13(e,g,:,:)=co';
                    val13(e,g,:,:)=-obj.kappa.gauss.weight(g)*bchi.M{e,g}'*dgamma0*obj.kappa.M{e,g}*obj.GammaC.mass{end}(e);

                    [ro,co]=meshgrid(obj.kappa.ind(e),bchi.ind(e));
                    row31(e,g,:,:)=ro';
                    col31(e,g,:,:)=co'; 
                    val31(e,g,:,:)=obj.kappa.gauss.weight(g)*obj.dt*obj.kappa.M{e,g}'*dgamma0'*bchi.M{e,g}*obj.GammaC.mass{end}(e);
                end
            end
            out13=sparse(row13(:),col13(:),val13(:));
            out13=bchi.tr'*out13;
            gammadof=accumarray(row3(:),val3(:));
            gammadof1=sparse(row31(:),col31(:),val31(:));
            gammadof1=gammadof1*bchi.tr;

            out31=obj.dxFB(gammadof,obj.kappa.dof)*gammadof1;
            out33=obj.dyFB(gammadof,obj.kappa.dof);

            out22=sparse(obj.r.Ndof,obj.r.Ndof);
            out23=sparse(obj.r.Ndof,obj.kappa.Ndof);
            out32=sparse(obj.kappa.Ndof,obj.r.Ndof);

            out11=out11/obj.unit.force*obj.unit.velocity;
            out12=out12/obj.unit.force*obj.unit.pressure;
            out13=out13/obj.unit.force*obj.unit.pressure;
            out21=out21/obj.unit.force*obj.unit.velocity;
            out22=out22/obj.unit.force*obj.unit.pressure;
            out23=out23/obj.unit.force*obj.unit.pressure;
            out31=out31*obj.unit.velocity;
            out32=out32*obj.unit.pressure;
            out33=out33*obj.unit.pressure;

            out=[out11,out12,out13;out21,out22,out23;out31,out32,out33];
        end

        function solve(obj)
            obj.init;

            fprintf("\nstart analysis\n");

            obj.init();
            good=0;

            while obj.t<obj.T
                % obj.v.dof=0;
                v0=obj.v.dof;
                chi0=obj.chi.dof;
                r0=obj.r.dof;
                kappa0=obj.kappa.dof;
                t0=obj.t;

                if obj.t+obj.dt>obj.T
                    obj.setting.step.dt=1.001*(obj.T-obj.t);
                end

                obj.t=obj.t+obj.dt;
                iter=0;
                retry=false;
                obj.computegap;
                res=obj.residual;
                stiff=obj.stiffness;

                fprintf('\nstep %d at time %.2f/%.2f with time increment %.2e\n',obj.setting.step.step+1,obj.t,obj.T,obj.dt);
                fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,norm(res),condest(stiff));

                while norm(res)>obj.setting.newton.tol
                    
                    if iter==obj.setting.newton.maxiter || condest(stiff)>1e12
                        retry=true;
                        break;
                    end

                    v1=obj.v.dof;
                    r1=obj.r.dof;
                    kappa1=obj.kappa.dof;

                    u=-stiff\res;

                    dv=obj.unit.velocity*u(1:obj.v.Ndof);
                    dr=obj.unit.pressure*u(obj.v.Ndof+1:obj.v.Ndof+obj.r.Ndof);
                    dkappa=obj.unit.pressure*u(end-obj.kappa.Ndof+1:end);
                    
                    eta=1;

                    obj.v.dof=v1+eta*dv;
                    obj.chi.dof=chi0+obj.dt*obj.v.dof;
                    obj.r.dof=r1+eta*dr;
                    obj.kappa.dof=kappa1+eta*dkappa;

                    while norm(obj.residual)>norm(res) && eta>obj.setting.newton.mineta
                        eta=eta/1.5;

                        obj.v.dof=v1+eta*dv;
                        obj.chi.dof=chi0+obj.dt*obj.v.dof;
                        obj.r.dof=r1+eta*dr;
                        obj.kappa.dof=kappa1+eta*dkappa;

                        fprintf('#');
                    end

                    obj.computegap;
                    res=obj.residual;
                    stiff=obj.stiffness;

                    iter=iter+1;

                    fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,norm(res),condest(stiff));
                end

                if norm(obj.dt*(obj.v.dof-v0),Inf)/obj.unit.velocity>obj.setting.tol.chi
                    retry=true;
                end

                if norm(obj.r.dof-r0,Inf)/obj.unit.pressure>obj.setting.tol.r
                    retry=true;
                end

                if norm(obj.kappa.dof-kappa0,Inf)/obj.unit.pressure>obj.setting.tol.kappa
                    retry=true;
                end

                if retry==true
                    obj.v.dof=v0;
                    obj.chi.dof=chi0;
                    obj.r.dof=r0;
                    obj.kappa.dof=kappa0;
                    obj.t=t0;
                  
                    obj.setting.step.dt=Inf;

                    if obj.setting.step.lev-1<-obj.setting.step.nmin
                        fprintf('\nanalysis failed\n');
                        return;
                    end

                    obj.setting.step.lev=obj.setting.step.lev-1;
                    fprintf("step failed\n");

                    good=0;

                    continue;
                end

                obj.setting.step.step=obj.setting.step.step+1;
                good=good+1;

                obj.write();

                fprintf("step completed\n");

                if good>=obj.setting.step.maxgood && obj.setting.step.lev+1<=obj.setting.step.nmax
                    obj.setting.step.lev=obj.setting.step.lev+1;
                end
            end

            fprintf("\nanalysis completed\n");
        end

        % other functions

        function out=dt(obj)
            if obj.setting.step.lev>=0
                p=(obj.setting.step.dtmax/obj.setting.step.dt0)^(1/obj.setting.step.nmax);
            else
                p=(obj.setting.step.dt0/obj.setting.step.dtmin)^(1/obj.setting.step.nmin);
            end

            out=p^obj.setting.step.lev*obj.setting.step.dt0;

            if obj.setting.step.dt<Inf
                out=obj.setting.step.dt;
                return;
            end
        end

        function computegap(obj)
            for a=1:obj.GammaMC.num(1)
                obj.GammaMC.dcoord(a,:)=obj.chiM(obj.GammaMC.coord(a,:)',obj.t);
                obj.GammaMC.compute;
            end

            for a=1:obj.OmegaM.num(1)
                obj.OmegaM.dcoord(a,:)=obj.chiM(obj.OmegaM.coord(a,:)',obj.t);
                %obj.OmegaM.compute;
            end

            bchi=obj.chi.trace(obj.GammaC,obj.r.gauss);

            obj.gamma=zeros(obj.GammaC.num(end),obj.kappa.gauss.num);
            obj.dgamma=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,2);

            for e=1:obj.GammaC.num(end)
                for gauss=1:obj.r.gauss.num
                    x=bchi.eval(e,gauss);

                    in=obj.OmegaM.in(x);
                    dm=Inf;

                    for a=1:obj.GammaMC.num(end)
                        x0=obj.GammaMC.dcoord(obj.GammaMC.topol{end}(a,:),:);
                        x1=x0(2,:)';
                        x0=x0(1,:)';
                        dx=x1-x0;
                        l=norm(dx);
                        n=obj.GammaMC.normal(a,:)';

                        alpha=zeros(1,3);
                        alpha(2)=1;
                        alpha(3)=(x-x0)'*dx/l^2;
                        if alpha(3)<0 || alpha(3)>1
                            alpha(3)=0;
                        end

                        x2=(1-alpha).*x0+alpha.*x1;
                        g=-x2+x;

                        d=min(vecnorm(g));

                        if d<dm
                            dm=d;
                            obj.gamma(e,gauss)=dm*(-in+~in);
                            obj.dgamma(e,gauss,:)=(eye(2)-1/l^2*(dx*dx'))*n;
                        end
                    end
                end   
            end
        end

        function computegapold(obj)
            for a=1:obj.GammaMC.num(1)
                obj.GammaMC.dcoord(a,:)=obj.chiM(obj.GammaMC.coord(a,:)',obj.t);
                obj.GammaMC.compute;
            end

            bchi=obj.chi.trace(obj.GammaC,obj.r.gauss);

            obj.gamma=zeros(obj.GammaC.num(end),obj.kappa.gauss.num);
            obj.dgamma=zeros(obj.GammaC.num(end),obj.kappa.gauss.num,2);

            for e=1:obj.GammaC.num(end)
                for gauss=1:obj.r.gauss.num
                    x=bchi.eval(e,gauss);

                    found=0;

                    for a=1:obj.GammaMC.num(end)
                        x0=obj.GammaMC.dcoord(obj.GammaMC.topol{end}(a,:),:);
                        x1=x0(2,:)';
                        x0=x0(1,:)';
                        dx=x1-x0;
                        l=norm(dx);
                        n=obj.GammaMC.normal(a,:)';

                        alpha=(x-x0)'*dx/l^2;

                        if alpha>=-1e-10 && alpha<=1+1e-10
                            obj.gamma(e,gauss)=(x-x0-alpha*dx)'*n;
                            obj.dgamma(e,gauss,:)=(eye(2)-1/l^2*(dx*dx'))*n;
                            found=found+1;
                        end
                    end

                    assert(found>0);
                end   
            end
        end

        % Ficher-Burmistaing

        function out=FB(obj,x,y)
            out=obj.setting.FB.cx*x+obj.setting.FB.cy*y-obj.ab(obj.setting.FB.cx*x-obj.setting.FB.cy*y);
        end

        function out=dxFB(obj,x,y)
            out=obj.setting.FB.cx*(1-obj.dab(obj.setting.FB.cx*x-obj.setting.FB.cy*y));
            n=length(out);
            out=sparse(1:n,1:n,out);
        end

        function out=dyFB(obj,x,y)
            out=obj.setting.FB.cy*(1+obj.dab(obj.setting.FB.cx*x-obj.setting.FB.cy*y));
            n=length(out);
            out=sparse(1:n,1:n,out);
        end

        function out=ab(obj,t)
            out=abs(t).^(1+obj.setting.FB.p);
        end

        function out=dab(obj,t)
            out=(1+obj.setting.FB.p)*sign(t).*abs(t).^obj.setting.FB.p;
        end
    end
end