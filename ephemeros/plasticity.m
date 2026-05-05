classdef plasticity < handle
    properties
        geo
        Omega
        GammaD
        GammaDn
        Gamma

        spc
        mlt

        unit

        unk
        mat
        dat

        setting

        t
        T
    end
    methods
        % initialize, save in files and manage external actions

        function init(obj)
            obj.spc.U=space(obj.Omega,2,type="Vec");
            obj.spc.R=space(obj.GammaD,0,type="Vec");
            obj.spc.Rn=space(obj.GammaDn,0);
            obj.spc.M=space(obj.Omega.skel,0,type="Vec");
            obj.spc.Z=space(obj.Omega,1,type="DevSym");
            obj.spc.L=space(obj.Omega,0);
            obj.spc.P=space(obj.Omega,0);

            obj.unk.u=field(obj.spc.U);
            obj.unk.dotu=field(obj.spc.U);
            obj.unk.r=field(obj.spc.R);
            obj.unk.rn=field(obj.spc.Rn);
            obj.unk.m=field(obj.spc.M);
            obj.unk.z=field(obj.spc.Z);
            obj.unk.dotz=field(obj.spc.Z);
            obj.unk.la=field(obj.spc.L);
            obj.unk.p=field(obj.spc.P);

            obj.t=0;
            obj.setting.step.dt=Inf;
            obj.setting.step.lev=0;
            obj.setting.step.step=0;

            obj.unit.time=obj.T;
            obj.unit.length=obj.geo.diamest;
            obj.unit.pressure=obj.mat.C.mu*obj.unit.strain;
            
            obj.unit.velocity=obj.unit.length/obj.unit.time;

            obj.mlt=multi;
            obj.mlt.add(obj.spc.U,scale=[obj.unit.pressure*obj.unit.length^(obj.geo.d-1),obj.unit.strain*obj.unit.length/obj.unit.time]);
            obj.mlt.add(obj.spc.R,scale=[obj.unit.length^obj.geo.d,obj.unit.pressure]);
            obj.mlt.add(obj.spc.M,scale=[obj.unit.length^(obj.geo.d-1),obj.unit.pressure*obj.unit.length]);
            obj.mlt.add(obj.spc.Z,scale=[obj.unit.pressure*obj.unit.length^obj.geo.d,obj.unit.strain/obj.unit.time]);
            obj.mlt.add(obj.spc.L,scale=[obj.unit.pressure*obj.unit.length^obj.geo.d,obj.unit.strain]);
            obj.mlt.add(obj.spc.P,scale=[obj.unit.strain*obj.unit.length^obj.geo.d,obj.unit.pressure]);
            obj.mlt.add(obj.spc.Rn,scale=[obj.unit.length^obj.geo.d,obj.unit.pressure]);
            
            I=tensorprod(eye(obj.geo.d),eye(obj.geo.d));
            Id=reshape(eye(obj.geo.d^2),[obj.geo.d,obj.geo.d,obj.geo.d,obj.geo.d]);
            Tr=permute(Id,[1,4,2,3]);
            obj.mat.C.ten=obj.mat.C.la*I+obj.mat.C.mu*(Id+Tr);
            obj.mat.D.ten=obj.mat.D.la*I+obj.mat.D.mu*(Id+Tr);
        end

        function write(obj)
            figure(1);
            clf;
            theme(gcf,"light");
            % obj.unk.u.plot(prc=@(u) u(1,:),trisurf=0,LineStyle="None");
            % colorbar;
            % subplot(1,2,2);
            obj.unk.z.plot(prc="norm",trisurf=0,LineStyle="None",displ=obj.unk.u,scale=100);
            colorbar;
            title(sprintf("t=%.4f",obj.t));
            drawnow;
            %pause; 

            % figure(2);
            % clf;
            % z0=obj.unk.z.dof;
            % obj.unk.z.dof=1/obj.unit.strain*obj.unk.z.dof;
            % obj.unk.z.plot(prc="norm",trisurf=1,LineStyle="None");
            % obj.unk.z.dof=z0;
            % colorbar;
            % title(sprintf("t=%.4f",obj.t));
            % drawnow;
        end

        % residual, stiffness and solver

        function [res,stiff]=assemble(obj)
            d=obj.geo.d;
            n=obj.geo.gauss{2}.trace{1}.normal;

            U=obj.spc.U;
            R=obj.spc.R;
            Rn=obj.spc.Rn;
            M=obj.spc.M;
            Z=obj.spc.Z;
            L=obj.spc.L;
            P=obj.spc.P;

            nU=U.ref.ndof;
            nR=R.ref.ndof;
            nRn=Rn.ref.ndof;
            nM=M.ref.ndof;
            nZ=Z.ref.ndof;
            nL=L.ref.ndof;
            nP=P.ref.ndof;

            sU=U.shape;
            sR=R.shape;
            sRn=Rn.shape;
            sM=M.shape;
            sZ=Z.shape;
            sL=L.shape;
            sP=P.shape;

            u=obj.unk.u;
            dotu=obj.unk.dotu;
            r=obj.unk.r;
            rn=obj.unk.rn;
            m=obj.unk.m;
            z=obj.unk.z;
            dotz=obj.unk.dotz;
            la=obj.unk.la;
            p=obj.unk.p;

            u.eval;
            dotu.eval;
            r.eval;
            rn.eval;
            m.eval;
            z.eval;
            dotz.eval;
            la.eval;
            p.eval;

            C=obj.mat.C.ten;
            D=obj.mat.D.ten;
            eta=obj.mat.eta;
            l=obj.mat.l;
            A=obj.mat.A;
            B=obj.mat.B;

            f=field(U);
            f.map(@(x) obj.dat.f(x,obj.t));
            f.eval;
            f=f.val{2,1};
            g=field(space(obj.Omega.bound,0,type="Vec"));
            g.map(@(x) obj.dat.g(x,obj.t));
            g.eval;

            TT=eta*tensorprod(D,dotu.val{2,2},[3,4],[1,2]);
            K=l^2*A*(eta*dotu.val{2,3}+u.val{2,3});

            Te=tensorprod(C,u.val{2,2},[3,4],[1,2]);
            Te=Te-tensorprod(C,z.val{2,1},[3,4],[1,2]);
            Te=Te-squeeze(tensorprod(tensorprod(C,1/d*eye(d),[3,4],[1,2]),la.val{2,1}));

            [Y,dY]=obj.mat.yield(p.val{2,1});
            [dG,d2G]=obj.mat.diss(dotz.val{2,1});

            Tp=B*(eta*dotz.val{2,1}+z.val{2,1})+reshape(Y,1,1,[]).*dG;
            Kp=l^2*B*(eta*dotz.val{2,2}+z.val{2,2});
            
            int=squeeze(sum(reshape(Te+TT,d,d,1,[]).*reshape(sU{2,2},d,d,nU,[]),[1,2]));
            int=int+squeeze(sum(reshape(K,d,d,d,1,[]).*reshape(sU{2,3},d,d,d,nU,[]),[1,2,3]));
            int=int-squeeze(sum(reshape(f,d,1,[]).*reshape(sU{2,1},d,nU,[]),1));
            res1=asmbvec(int,2,2,2,U);
            int=-squeeze(sum(reshape(r.val{2,1},d,1,[]).*reshape(sU{1,1},d,nU,[]),1));
            res1=res1+asmbvec(int,2,1,1,U);
            int=-squeeze(sum(reshape(n.*reshape(rn.val{2,1},1,[]),d,1,[]).*sU{1,1},1));
            res1=res1+asmbvec(int,2,1,1,U);
            int=squeeze(sum(squeeze(sum(reshape(sU{1,2},d,d,nU,[]).*reshape(n,1,d,1,[]),2)).*reshape(m.val{2,1},d,1,[]),1));
            res1=res1+asmbvec(int,2,1,1,U);

            int=-squeeze(sum(reshape(u.val{1,1}-g.val{2,1},d,1,[]).*reshape(sR{2,1},d,nR,[]),1));
            res2=asmbvec(int,2,1,2,R);

            int=-squeeze(sum(reshape(sum((u.val{1,1}-g.val{2,1}).*n,1),1,[]).*reshape(sRn{2,1},nRn,[]),1));
            res7=asmbvec(int,2,1,2,Rn);

            int=squeeze(sum(reshape(squeeze(sum(reshape(u.val{1,2},d,d,[]).*reshape(n,1,d,[]),2)),d,1,[]).*reshape(sM{2,1},d,nM,[]),1));
            res3=asmbvec(int,2,1,2,M);

            int=squeeze(sum(reshape(Tp-Te,d,d,1,[]).*reshape(sZ{2,1},d,d,nZ,[]),[1,2]));
            int=int+squeeze(sum(reshape(Kp,d,d,d,1,[]).*reshape(sZ{2,2},d,d,d,nZ,[]),[1,2,3]));
            res4=asmbvec(int,2,2,2,Z);

            int=-squeeze(sum(reshape(tensorprod(1/d*eye(d),Te,[1,2],[1,2]),1,[]).*reshape(sL{2,1},nL,[]),1));
            int=int-squeeze(sum(reshape(p.val{2,1},1,[]).*reshape(sL{2,1},nL,[]),1));
            res5=asmbvec(int,2,2,2,L);

            int=-squeeze(sum(reshape(la.val{2,1},1,[]).*reshape(sP{2,1},nP,[]),1));
            res6=asmbvec(int,2,2,2,L);

            obj.mlt.rescomp={res1,res2,res3,res4,res5,res6,res7};

            int=squeeze(sum(reshape(tensorprod(eta*D+obj.dt*C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(sU{2,2},d,d,1,nU,[]),[1,2]));
            int=int+l^2*A*(eta+obj.dt)*squeeze(sum(reshape(sU{2,3},d,d,d,nU,1,[]).*reshape(sU{2,3},d,d,d,1,nU,[]),[1,2,3]));
            stiff11=asmbmat(int,2,2,2,2,U,U);

            int=-squeeze(sum(reshape(sU{1,1},d,nU,1,[]).*reshape(sR{2,1},d,1,nR,[]),1));
            stiff12=asmbmat(int,2,1,1,2,U,R);

            int=-reshape(reshape(sum(reshape(n,d,1,[]).*sU{1,1},1),nU,1,[]).*reshape(sRn{2,1},1,nRn,[]),nU,nRn,[]);
            stiff17=asmbmat(int,2,1,1,2,U,Rn);

            int=squeeze(sum(reshape(squeeze(sum(sU{1,2}.*reshape(n,1,d,1,[]),2)),d,nU,1,[]).*reshape(sM{2,1},d,1,nM,[]),1));
            stiff13=asmbmat(int,2,1,1,2,U,M);

            int=-obj.dt*squeeze(sum(reshape(tensorprod(C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(sZ{2,1},d,d,1,nZ,[]),[1,2]));
            stiff14=asmbmat(int,2,2,2,2,U,Z);

            int=-sum(reshape(tensorprod(C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(tensorprod(1/d*eye(d),sL{2,1}),d,d,1,nL,[]),[1,2]);
            stiff15=asmbmat(int,2,2,2,2,U,L);

            int=obj.dt*sum(reshape(tensorprod(C,sZ{2,1},[3,4],[1,2]),d,d,nZ,1,[]).*reshape(sZ{2,1},d,d,1,nZ,[]),[1,2]);
            int=int+B*(eta+obj.dt)*sum(reshape(sZ{2,1},d,d,nZ,1,[]).*reshape(sZ{2,1},d,d,1,nZ,[]),[1,2]);
            int=squeeze(int)+squeeze(sum(reshape(squeeze(sum(reshape((reshape(Y,1,1,1,1,[]).*d2G),d,d,d,d,1,[]).*reshape(sZ{2,1},1,1,d,d,nZ,[]),[3,4])),d,d,nZ,1,[]).*reshape(sZ{2,1},d,d,1,nZ,[]),[1,2]));
            int=int+squeeze(l^2*B*(eta+obj.dt)*sum(reshape(sZ{2,2},d,d,d,nZ,1,[]).*reshape(sZ{2,2},d,d,d,1,nZ,[]),[1,2,3]));
            stiff44=asmbmat(int,2,2,2,2,Z,Z);

            int=squeeze(sum(reshape(sZ{2,1},d,d,nZ,1,[]).*reshape(tensorprod(tensorprod(C,1/d*eye(d),[3,4],[1,2]),sL{2,1}),d,d,1,nL,[]),[1,2]));
            stiff45=asmbmat(int,2,2,2,2,Z,L);

            int=squeeze(reshape(sum(reshape(reshape(dY,1,1,[]).*dG,d,d,1,[]).*sZ{2,1},[1,2]),nZ,1,[]).*reshape(sP{2,1},1,nP,[]));
            stiff46=asmbmat(int,2,2,2,2,Z,P);

            int=1/d^2*tensorprod(tensorprod(C,eye(d),[3,4],[1,2]),eye(d),[1,2],[1,2])*reshape(sL{2,1},nL,1,[]).*reshape(sL{2,1},1,nL,[]);
            stiff55=asmbmat(int,2,2,2,2,L,L);

            int=-reshape(sL{2,1},nL,1,[]).*reshape(sP{2,1},1,nL,[]);
            stiff56=asmbmat(int,2,2,2,2,L,P);

            obj.mlt.stiffblk={stiff11,stiff12,stiff13,stiff14,stiff15,[],stiff17;
                obj.dt*stiff12',[],[],[],[],[],[];
                obj.dt*stiff13',[],[],[],[],[],[];
                stiff14',[],[],stiff44,stiff45,stiff46,[];
                obj.dt*stiff15',[],[],obj.dt*stiff45',stiff55,stiff56,[];
                [],[],[],[],stiff56',[],[];
                obj.dt*stiff17',[],[],[],[],[],[]};

            res=obj.mlt.res;
            stiff=obj.mlt.stiff;
        end

        function solve(obj)
            obj.init;

            fprintf("\nstart analysis\n");

            obj.init();
            good=0;

            while obj.t<obj.T
                dotu0=obj.unk.dotu.dof;
                u0=obj.unk.u.dof;
                r0=obj.unk.r.dof;
                rn0=obj.unk.rn.dof;
                m0=obj.unk.m.dof;
                dotz0=obj.unk.dotz.dof;
                z0=obj.unk.z.dof;
                la0=obj.unk.la.dof;
                p0=obj.unk.p.dof;
                t0=obj.t;

                if obj.t+obj.dt>obj.T
                    obj.setting.step.dt=1.001*(obj.T-obj.t);
                end

                obj.t=obj.t+obj.dt;
                iter=0;
                retry=false;

                obj.unk.u.dof=u0+obj.dt*dotu0;
                obj.unk.z.dof=z0+obj.dt*dotz0;

                [sres,sstiff]=obj.assemble;
                %[P,R,C]=equilibrate(stiff);
                %sstiff=R*P*stiff*C;
                %sres=R*P*res;
                %con=condest(sstiff);
                con=-1;
                nor=norm(sres,Inf);
                nor0=nor;

                tdelta=zeros(size(sres));

                fprintf('\nstep %d at time %.4f/%.4f with time increment %.2e\n',obj.setting.step.step+1,obj.t,obj.T,obj.dt);
                fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/min([1,nor0])/obj.setting.newton.tol,con);

                while nor>min([1,nor0])*obj.setting.newton.tol
                    if iter==obj.setting.newton.maxiter  || con>1e12
                        retry=true;
                    end

                    tic

                    delta=-sstiff\sres;
                    %delta=C*delta;

                    toc

                    tic

                    ddotu=obj.mlt.comp(delta,1);
                    dr=obj.mlt.comp(delta,2);
                    drn=obj.mlt.comp(delta,7);
                    dm=obj.mlt.comp(delta,3);
                    ddotz=obj.mlt.comp(delta,4);
                    dla=obj.mlt.comp(delta,5);
                    dp=obj.mlt.comp(delta,6);

                    eta=1.5;

                    nor1=nor;
                    dotu1=obj.unk.dotu.dof;
                    r1=obj.unk.r.dof;
                    rn1=obj.unk.rn.dof;
                    m1=obj.unk.m.dof;
                    dotz1=obj.unk.dotz.dof;
                    la1=obj.unk.la.dof;
                    p1=obj.unk.p.dof;

                    while nor>=0.9*nor1 || eta==1.5
                        eta=eta/1.5;

                        obj.unk.dotu.dof=dotu1+eta*ddotu;
                        obj.unk.u.dof=u0+obj.dt*obj.unk.dotu.dof;
                        obj.unk.r.dof=r1+eta*dr;
                        obj.unk.rn.dof=rn1+eta*drn;
                        obj.unk.m.dof=m1+eta*dm;
                        obj.unk.dotz.dof=dotz1+eta*ddotz;
                        obj.unk.z.dof=z0+obj.dt*obj.unk.dotz.dof;
                        obj.unk.la.dof=la1+eta*dla;
                        obj.unk.p.dof=p1+eta*dp;

                        [sres,sstiff]=obj.assemble;
                        % sstiff=R*P*stiff*C;
                        % sres=R*P*res;
                        %con=condest(sstiff);
                        con=-1;
                        nor=norm(sres,Inf);

                        fprintf("#");

                        if eta<0.42
                            retry=true;
                            break;
                        end
                    end

                    tdelta=tdelta+delta*eta;

                    toc

                    iter=iter+1;

                    fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/min([1,nor0])/obj.setting.newton.tol,con);
                    %fprintf('eig=%.5e\n',eigs(sstiff,1,'smallestabs'));

                    if retry==true
                        break;
                    end
                end

                if obj.dt/obj.T*norm(tdelta,Inf)/max(1,norm(delta))>obj.setting.tol
                    retry=true;
                    fprintf("*");
                end

                if retry==true
                    obj.unk.v.dof=dotu0;
                    obj.unk.u.dof=u0;
                    obj.unk.r.dof=r0;
                    obj.unk.rn.dof=rn0;
                    obj.unk.m.dof=m0;
                    obj.unk.dotz.dof=dotz0;
                    obj.unk.z.dof=z0;
                    obj.unk.la.dof=la0;
                    obj.unk.p.dof=p0;
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
    end
end