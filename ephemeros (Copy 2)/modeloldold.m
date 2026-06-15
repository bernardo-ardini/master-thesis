classdef modeloldold < handle
    properties
        geo
        Omega
        GammaD
        GammaDn
        Gamma

        name

        spc
        mlt

        unit

        hist
        unk
        var
        mat
        dat

        tols
        toln
        told
        tolal
        tolpi

        step
        dt
        dtmax
        dt0
        dts
        dtd
        dtn
        t
        T

        scale
    end
    methods
        % initialize, save in files and manage external actions

        function init(obj)
            obj.spc.U=space(obj.Omega,2,type="Vec");
            obj.spc.R=space(obj.GammaD,0,type="Vec");
            obj.spc.Rn=space(obj.GammaDn,0);
            obj.spc.M=space(obj.Omega.skel,0,type="Vec");
            obj.spc.Pi=space(obj.Omega,1,type="DevSym");
            obj.spc.Al=space(obj.Omega,1);
            obj.spc.La=space(obj.Omega,0);
            obj.spc.Las=space(obj.Omega,0);

            obj.unk.u=field(obj.spc.U);
            obj.unk.dotu=field(obj.spc.U);
            obj.unk.r=field(obj.spc.R);
            obj.unk.rn=field(obj.spc.Rn);
            obj.unk.m=field(obj.spc.M);
            obj.unk.pi=field(obj.spc.Pi);
            obj.unk.dotpi=field(obj.spc.Pi);
            obj.unk.al=field(obj.spc.Al);
            obj.unk.dotal=field(obj.spc.Al);
            obj.unk.la=field(obj.spc.La);
            obj.unk.p=field(obj.spc.Las);

            obj.mlt=multi;
            obj.mlt.add(obj.spc.U);
            obj.mlt.add(obj.spc.R);
            obj.mlt.add(obj.spc.Rn);
            obj.mlt.add(obj.spc.M);
            obj.mlt.add(obj.spc.Pi);
            obj.mlt.add(obj.spc.Al);
            obj.mlt.add(obj.spc.La);
            obj.mlt.add(obj.spc.Las);

            figure(2);
            set(gca,"NextPlot","replacechildren");
            theme(gcf,"light");

            obj.hist={};
        end

        function plot(obj)
            clf;

            subplot(2,2,1);
            obj.unk.pi.plot(1,@(pi) pagenorm(pi,"fro"));
            colormap parula(30);
            colorbar;
            %clim([0,1.5e-3]);
            title("\pi");

            subplot(2,2,2);
            obj.unk.p.plot(1,@(p) p);
            colormap parula(30);
            colorbar;
            title("p");

            subplot(2,2,3);
            obj.unk.al.plot(1,@(al) al);
            colormap parula(30);
            colorbar;
            %clim([0,1.5e-3]);
            title("\alpha");

            subplot(2,2,4);
            obj.unk.la.plot(1,@(al) al,disp=obj.hist{end}.displ,scale=obj.scale);
            colormap parula(30);
            colorbar;
            %clim([0,1.5e-3]);
            title("\lambda");

            sgtitle(sprintf("t=%.4f",obj.t));
            drawnow;
        end

        % residual, stiffness and solver

        function [res,stiff,diss]=assemble(obj)
            d=obj.geo.d;
            n=obj.geo.gauss{2}.trace{1}.normal;

            U=obj.spc.U;
            R=obj.spc.R;
            Rn=obj.spc.Rn;
            M=obj.spc.M;
            Pi=obj.spc.Pi;
            Al=obj.spc.Al;
            La=obj.spc.La;
            Las=obj.spc.Las;

            nU=U.ref.ndof;
            nR=R.ref.ndof;
            nRn=Rn.ref.ndof;
            nM=M.ref.ndof;
            nPi=Pi.ref.ndof;
            nAl=Al.ref.ndof;
            nLa=La.ref.ndof;
            nLas=Las.ref.ndof;

            sU=U.shape;
            sR=R.shape;
            sRn=Rn.shape;
            sM=M.shape;
            sPi=Pi.shape;
            sAl=Al.shape;
            sLa=La.shape;
            sLas=Las.shape;

            u=obj.unk.u;
            dotu=obj.unk.dotu;
            r=obj.unk.r;
            rn=obj.unk.rn;
            m=obj.unk.m;
            pi=obj.unk.pi;
            dotpi=obj.unk.dotpi;
            al=obj.unk.al;
            dotal=obj.unk.dotal;
            la=obj.unk.la;
            p=obj.unk.p;

            u.eval;
            dotu.eval;
            r.eval;
            rn.eval;
            m.eval;
            pi.eval;
            dotpi.eval;
            al.eval;
            dotal.eval;
            la.eval;
            p.eval;

            a=obj.mat.a;
            b=obj.mat.b;
            c=obj.mat.c;
            C=obj.mat.C;

            D=obj.mat.D;
            dd=obj.mat.d;
            e=obj.mat.e;

            f=field(U);
            f.map(@(x) obj.dat.f(x,obj.t));
            f.eval;
            f=f.val{2,1};
            g=field(space(obj.Omega.bound,0,type="Vec"));
            g.map(@(x) obj.dat.g(x,obj.t));
            g.eval;

            TT=tensorprod(D,dotu.val{2,2},[3,4],[1,2]);
            obj.var.T=TT;
            K=a*u.val{2,3};
            Ep=pi.val{2,1}+1/d*reshape(tensorprod(eye(d),la.val{2,1}),d,d,[]);
            E=1/2*(u.val{2,2}+permute(u.val{2,2},[2,1,3]));
            Ee=E-Ep;
            Te=reshape(sum(reshape(C,d,d,d,d,1).*reshape(Ee,1,1,d,d,[]),[3,4]),d,d,[]);
            obj.var.Te=Te;
            [Y,dY]=obj.mat.Y(p.val{2,1});
            [dRR,d2RR]=obj.mat.R(dotpi.val{2,1});
            [drhoQ,dbeQ,drhorhoQ,drhobeQ,dbebeQ]=obj.mat.Q(dotpi.val{2,1},dotal.val{2,1});
            Tp=dd*dotpi.val{2,1}+drhoQ+reshape(Y,1,1,[]).*dRR;
            Kp=b*pi.val{2,2};
            [dWd,d2Wd]=obj.mat.Wd(al.val{2,1});
            tau=dWd+dbeQ+e*dotal.val{2,1};
            kappa=c*al.val{2,2};
            [L,dL,d2L]=obj.mat.L(al.val{2,1});
            
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
            res3=asmbvec(int,2,1,2,Rn);

            int=squeeze(sum(reshape(squeeze(sum(reshape(u.val{1,2},d,d,[]).*reshape(n,1,d,[]),2)),d,1,[]).*reshape(sM{2,1},d,nM,[]),1));
            res4=asmbvec(int,2,1,2,M);

            int=squeeze(sum(reshape(Tp-Te,d,d,1,[]).*reshape(sPi{2,1},d,d,nPi,[]),[1,2]));
            int=int+squeeze(sum(reshape(Kp,d,d,d,1,[]).*reshape(sPi{2,2},d,d,d,nPi,[]),[1,2,3]));
            res5=asmbvec(int,2,2,2,Pi);

            int=squeeze(reshape(dL.*p.val{2,1},1,[]).*reshape(sAl{2,1},nAl,[]));
            int=int+squeeze(reshape(tau,1,[]).*reshape(sAl{2,1},nAl,[]));
            int=int+squeeze(sum(reshape(kappa,d,1,[]).*reshape(sAl{2,2},d,nAl,[]),1));
            res6=asmbvec(int,2,2,2,Al);

            int=-squeeze(sum(reshape(tensorprod(1/d*eye(d),Te,[1,2],[1,2]),1,[]).*reshape(sLa{2,1},nLa,[]),1));
            int=int-squeeze(sum(reshape(p.val{2,1},1,[]).*reshape(sLa{2,1},nLa,[]),1));
            res7=asmbvec(int,2,2,2,La);

            int=-squeeze(sum(reshape(la.val{2,1},1,[]).*reshape(sLas{2,1},nLas,[]),1));
            int=int+squeeze(sum(reshape(L,1,[]).*reshape(sLas{2,1},nLas,[]),1));
            res8=asmbvec(int,2,2,2,Las);

            obj.mlt.rescomp={res1,res2,res3,res4,res5,res6,res7,res8};

            int=squeeze(sum(reshape(tensorprod(D+obj.dt*C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(sU{2,2},d,d,1,nU,[]),[1,2]));
            int=int+a*obj.dt*squeeze(sum(reshape(sU{2,3},d,d,d,nU,1,[]).*reshape(sU{2,3},d,d,d,1,nU,[]),[1,2,3]));
            stiff11=asmbmat(int,2,2,2,2,U,U);

            int=-squeeze(sum(reshape(sU{1,1},d,nU,1,[]).*reshape(sR{2,1},d,1,nR,[]),1));
            stiff12=asmbmat(int,2,1,1,2,U,R);

            int=-reshape(reshape(sum(reshape(n,d,1,[]).*sU{1,1},1),nU,1,[]).*reshape(sRn{2,1},1,nRn,[]),nU,nRn,[]);
            stiff13=asmbmat(int,2,1,1,2,U,Rn);

            int=squeeze(sum(reshape(squeeze(sum(sU{1,2}.*reshape(n,1,d,1,[]),2)),d,nU,1,[]).*reshape(sM{2,1},d,1,nM,[]),1));
            stiff14=asmbmat(int,2,1,1,2,U,M);

            int=-obj.dt*squeeze(sum(reshape(tensorprod(C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(sPi{2,1},d,d,1,nPi,[]),[1,2]));
            stiff15=asmbmat(int,2,2,2,2,U,Pi);

            int=-sum(reshape(tensorprod(C,sU{2,2},[3,4],[1,2]),d,d,nU,1,[]).*reshape(tensorprod(1/d*eye(d),sLa{2,1}),d,d,1,nLa,[]),[1,2]);
            stiff17=asmbmat(int,2,2,2,2,U,La);

            int=obj.dt*squeeze(sum(reshape(tensorprod(C,sPi{2,1},[3,4],[1,2]),d,d,nPi,1,[]).*reshape(sPi{2,1},d,d,1,nPi,[]),[1,2]));
            int=int+dd*squeeze(sum(reshape(sPi{2,1},d,d,nPi,1,[]).*reshape(sPi{2,1},d,d,1,nPi,[]),[1,2]));
            int=int+reshape(Y,1,1,[]).*squeeze(sum(reshape(sum(reshape(d2RR,d,d,d,d,1,[]).*reshape(sPi{2,1},1,1,d,d,nPi,[]),[3,4]),d,d,1,nPi,[]).*reshape(sPi{2,1},d,d,nPi,1,[]),[1,2]));
            int=int+squeeze(sum(reshape(sum(reshape(drhorhoQ,d,d,d,d,1,[]).*reshape(sPi{2,1},1,1,d,d,nPi,[]),[3,4]),d,d,1,nPi,[]).*reshape(sPi{2,1},d,d,nPi,1,[]),[1,2]));
            int=int+squeeze(b*obj.dt*sum(reshape(sPi{2,2},d,d,d,nPi,1,[]).*reshape(sPi{2,2},d,d,d,1,nPi,[]),[1,2,3]));
            stiff55=asmbmat(int,2,2,2,2,Pi,Pi);

            int=reshape(sum(reshape(drhobeQ,d,d,1,[]).*reshape(sPi{2,1},d,d,nPi,[]),[1,2]),nPi,1,[]).*reshape(sAl{2,1},1,nAl,[]);
            stiff56=asmbmat(int,2,2,2,2,Pi,Al);

            int=squeeze(sum(reshape(sPi{2,1},d,d,nPi,1,[]).*reshape(tensorprod(tensorprod(C,1/d*eye(d),[3,4],[1,2]),sLa{2,1}),d,d,1,nLa,[]),[1,2]));
            stiff57=asmbmat(int,2,2,2,2,Pi,La);

            int=reshape(dY,1,1,[]).*reshape(sum(reshape(dRR,d,d,1,[]).*sPi{2,1},[1,2]),nPi,1,[]).*reshape(sLas{2,1},1,nLas,[]);
            stiff58=asmbmat(int,2,2,2,2,Pi,Las);

            int=obj.dt*reshape(reshape(d2L,1,[]).*sAl{2,1},nAl,1,[]).*reshape(sAl{2,1},1,nAl,[]).*reshape(p.val{2,1},1,1,[]);
            int=int+obj.dt*reshape(d2Wd,1,1,[]).*reshape(sAl{2,1},nAl,1,[]).*reshape(sAl{2,1},1,nAl,[]);
            int=int+obj.mat.e*reshape(sAl{2,1},nAl,1,[]).*reshape(sAl{2,1},1,nAl,[]);
            int=int+obj.dt*obj.mat.c*squeeze(sum(reshape(sAl{2,2},d,nAl,1,[]).*reshape(sAl{2,2},d,1,nAl,[]),1));
            int=int+reshape(dbebeQ,1,1,[]).*reshape(sAl{2,1},nAl,1,[]).*reshape(sAl{2,1},1,nAl,[]);
            stiff66=asmbmat(int,2,2,2,2,Al,Al);

            int=reshape(dL,1,1,[]).*reshape(sAl{2,1},nAl,1,[]).*reshape(sLas{2,1},1,nLas,[]);
            stiff68=asmbmat(int,2,2,2,2,Al,Las);

            int=1/d^2*tensorprod(tensorprod(C,eye(d),[3,4],[1,2]),eye(d),[1,2],[1,2])*reshape(sLa{2,1},nLa,1,[]).*reshape(sLa{2,1},1,nLa,[]);
            stiff77=asmbmat(int,2,2,2,2,La,La);

            int=-reshape(sLa{2,1},nLa,1,[]).*reshape(sLas{2,1},1,nLa,[]);
            stiff78=asmbmat(int,2,2,2,2,La,Las);

            obj.mlt.stiffblk={stiff11,stiff12,stiff13,stiff14,stiff15,[],stiff17,[];
                obj.dt*stiff12',[],[],[],[],[],[],[];
                obj.dt*stiff13',[],[],[],[],[],[],[];
                obj.dt*stiff14',[],[],[],[],[],[],[];
                stiff15',[],[],[],stiff55,stiff56,stiff57,stiff58;
                [],[],[],[],stiff56',stiff66,[],stiff68;
                obj.dt*stiff17',[],[],[],obj.dt*stiff57',[],stiff77,stiff78;
                [],[],[],[],[],obj.dt*stiff68',stiff78,[]};

            res=obj.mlt.res;
            stiff=obj.mlt.stiff;

            diss=reshape(sum(reshape(sum(reshape(D,d,d,d,d,1).*reshape(dotu.val{2,2},1,1,d,d,[]),[3,4]),d,d,[]).*dotu.val{2,2},[1,2]),1,[])*reshape(obj.geo.gauss{2}.trace{2}.measure,[],1);
            diss=diss+dd*reshape(sum(dotpi.val{2,1}.*dotpi.val{2,1},[1,2]),1,[])*reshape(obj.geo.gauss{2}.trace{2}.measure,[],1);
            diss=diss+e*reshape(dotal.val{2,1}.*dotal.val{2,1},1,[])*reshape(obj.geo.gauss{2}.trace{2}.measure,[],1);
            diss=0.5*diss*obj.dt;
        end

        function solve(obj)
            fprintf("\nstart analysis\n");

            obj.init();
            
            obj.t=0;
            obj.step=0;
            obj.dt=obj.dt0;
            obj.dts=obj.dt0;
            obj.dtd=obj.dt0;
            obj.dtn=obj.dt0;

            while obj.t<obj.T
                dotu0=obj.unk.dotu.dof;
                u0=obj.unk.u.dof;
                r0=obj.unk.r.dof;
                rn0=obj.unk.rn.dof;
                m0=obj.unk.m.dof;
                dotpi0=obj.unk.dotpi.dof;
                pi0=obj.unk.pi.dof;
                dotal0=obj.unk.dotal.dof;
                al0=obj.unk.al.dof;
                la0=obj.unk.la.dof;
                p0=obj.unk.p.dof;
                t0=obj.t;

                obj.dt=min([obj.dts,obj.dtd,obj.dtn,1.5*obj.dt,obj.dtmax]);

                if obj.t+obj.dt>obj.T
                    obj.dt=(1+1e-5)*(obj.T-obj.t);
                end

                obj.t=obj.t+obj.dt;
                iter=0;
                retry=false;

                obj.unk.u.dof=u0+obj.dt*dotu0;
                obj.unk.pi.dof=pi0+obj.dt*dotpi0;
                obj.unk.al.dof=al0+obj.dt*dotal0;

                [sres,sstiff,diss]=obj.assemble;
                con=-1;
                nor=norm(sres)/sqrt(length(sres));

                tdel=0*[dotu0;dotpi0;dotal0];
                tdelpi=0*dotpi0;
                tdelal=0*dotal0;

                fprintf('\nstep %d at time %.4f/%.4f with time increment %.2e\n',obj.step+1,obj.t,obj.T,obj.dt);
                fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/obj.toln,con);

                while nor>obj.toln || iter==0
                    if iter==7
                        retry=true;
                        break;
                    end

                    delta=-sstiff\sres;
                    %delta=C*delta;

                    ddotu=obj.mlt.comp(delta,1);
                    dr=obj.mlt.comp(delta,2);
                    drn=obj.mlt.comp(delta,3);
                    dm=obj.mlt.comp(delta,4);
                    ddotpi=obj.mlt.comp(delta,5);
                    ddotal=obj.mlt.comp(delta,6);
                    dla=obj.mlt.comp(delta,7);
                    dp=obj.mlt.comp(delta,8);

                    tdel=tdel+[ddotu;ddotpi;ddotal];
                    tdelpi=tdelpi+ddotpi;
                    tdelal=tdelal+ddotal;

                    nor1=nor;
                    dotu1=obj.unk.dotu.dof;
                    r1=obj.unk.r.dof;
                    rn1=obj.unk.rn.dof;
                    m1=obj.unk.m.dof;
                    dotpi1=obj.unk.dotpi.dof;
                    dotal1=obj.unk.dotal.dof;
                    la1=obj.unk.la.dof;
                    p1=obj.unk.p.dof;

                    eta=1;
                    c=1e-4;

                    [sres,sstiff,diss]=obj.assemble;
                    con=-1;
                    nor=norm(sres)/sqrt(length(sres));

                    while nor>(1-eta*c)*nor1 && eta>8e-3
                        obj.unk.dotu.dof=dotu1+eta*ddotu;
                        obj.unk.u.dof=u0+obj.dt*obj.unk.dotu.dof;
                        obj.unk.r.dof=r1+eta*dr;
                        obj.unk.rn.dof=rn1+eta*drn;
                        obj.unk.m.dof=m1+eta*dm;
                        obj.unk.dotpi.dof=dotpi1+eta*ddotpi;
                        obj.unk.pi.dof=pi0+obj.dt*obj.unk.dotpi.dof;
                        obj.unk.dotal.dof=dotal1+eta*ddotal;
                        obj.unk.al.dof=al0+obj.dt*obj.unk.dotal.dof;
                        obj.unk.la.dof=la1+eta*dla;
                        obj.unk.p.dof=p1+eta*dp;

                        [sres,sstiff,diss]=obj.assemble;
                        con=-1;
                        nor=norm(sres)/sqrt(length(sres));

                        fprintf("*");
                        eta=eta/2;
                    end

                    iter=iter+1;

                    fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/obj.toln,con);
                end

                if obj.dt*norm(tdelal,Inf)>obj.tolal || obj.dt*norm(tdelpi,Inf)>obj.tolpi
                    retry=true;
                end

                if retry
                    obj.unk.dotu.dof=dotu0;
                    obj.unk.u.dof=u0;
                    obj.unk.r.dof=r0;
                    obj.unk.rn.dof=rn0;
                    obj.unk.m.dof=m0;
                    obj.unk.dotpi.dof=dotpi0;
                    obj.unk.pi.dof=pi0;
                    obj.unk.dotal.dof=dotal0;
                    obj.unk.al.dof=al0;
                    obj.unk.la.dof=la0;
                    obj.unk.p.dof=p0;
                    obj.t=t0;

                    obj.dtn=obj.dt/2;

                    fprintf("step failed\n");

                    continue;
                end

                err=norm(obj.dt*tdel)/sqrt(length(tdel));
                obj.dts=obj.dt*sqrt(obj.tols/err);

                if diss>1.1*obj.told
                    obj.dtd=obj.dt/1.5;
                else
                    obj.dtd=Inf;
                end          

                if iter>=4
                    obj.dtn=obj.dt/1.5;
                else
                    obj.dtn=Inf;
                end   

                obj.step=obj.step+1;

                h.t=obj.t;
                h.u=obj.unk.u.val;
                h.al=obj.unk.al.val;
                h.pi=obj.unk.pi.val;
                h.la=obj.unk.la.val;
                h.p=obj.unk.p.val;
                h.Te=obj.var.Te;
                h.t=obj.t;
                h.displ.coord=obj.unk.u.space.coord;
                h.displ.dof=obj.unk.u.dof;
                h.displ.num=obj.unk.u.space.num;
                h.displ.dim=obj.unk.u.space.basis.dim;

                obj.hist{end+1}=h;

                obj.plot();

                fprintf("step completed with dissipation %e\n",diss);

                %obj.setting.step.dt=min([obj.setting.step.dt*(obj.setting.diss/(diss+1e-9)),obj.setting.step.dtmax]);
            end

            fprintf("\nanalysis completed\n");
        end

        function oldsolve(obj)
            fprintf("\nstart analysis\n");

            obj.init();
            good=0;

            while obj.t<obj.T
                dotu0=obj.unk.dotu.dof;
                u0=obj.unk.u.dof;
                r0=obj.unk.r.dof;
                rn0=obj.unk.rn.dof;
                m0=obj.unk.m.dof;
                dotpi0=obj.unk.dotpi.dof;
                pi0=obj.unk.pi.dof;
                dotal0=obj.unk.dotal.dof;
                al0=obj.unk.al.dof;
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
                obj.unk.pi.dof=pi0+obj.dt*dotpi0;
                obj.unk.al.dof=al0+obj.dt*dotal0;

                [sres,sstiff,~]=obj.assemble;
                con=-1;
                nor=norm(sres)/sqrt(length(sres));
                nor0=nor;

                tdelta=zeros(size(sres));
                delta=zeros(size(sres));

                fprintf('\nstep %d at time %.4f/%.4f with time increment %.2e\n',obj.setting.step.step+1,obj.t,obj.T,obj.dt);
                fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/max([1,nor0])/obj.setting.newton.tol,con);

                while nor>max([1,nor0])*obj.setting.newton.tol
                    if iter==obj.setting.newton.maxiter || con>1e12
                        retry=true;
                        break;
                    end

                    delta=-sstiff\sres;
                    %delta=C*delta;

                    ddotu=obj.mlt.comp(delta,1);
                    dr=obj.mlt.comp(delta,2);
                    drn=obj.mlt.comp(delta,3);
                    dm=obj.mlt.comp(delta,4);
                    ddotpi=obj.mlt.comp(delta,5);
                    ddotal=obj.mlt.comp(delta,6);
                    dla=obj.mlt.comp(delta,7);
                    dp=obj.mlt.comp(delta,8);

                    nor1=nor;
                    dotu1=obj.unk.dotu.dof;
                    r1=obj.unk.r.dof;
                    rn1=obj.unk.rn.dof;
                    m1=obj.unk.m.dof;
                    dotpi1=obj.unk.dotpi.dof;
                    dotal1=obj.unk.dotal.dof;
                    la1=obj.unk.la.dof;
                    p1=obj.unk.p.dof;

                    tdelta=tdelta+delta;
                    eta=1;

                    while nor>(1-eta*c)*norm1
                        obj.unk.dotu.dof=dotu1+eta*ddotu;
                        obj.unk.u.dof=u0+obj.dt*obj.unk.dotu.dof;
                        obj.unk.r.dof=r1+eta*dr;
                        obj.unk.rn.dof=rn1+eta*drn;
                        obj.unk.m.dof=m1+eta*dm;
                        obj.unk.dotpi.dof=dotpi1+eta*ddotpi;
                        obj.unk.pi.dof=pi0+obj.dt*obj.unk.dotpi.dof;
                        obj.unk.dotal.dof=dotal1+eta*ddotal;
                        obj.unk.al.dof=al0+obj.dt*obj.unk.dotal.dof;
                        obj.unk.la.dof=la1+eta*dla;
                        obj.unk.p.dof=p1+eta*dp;
                    end                    

                    iter=iter+1;

                    [sres,sstiff,~]=obj.assemble;
                    con=-1;
                    nor=norm(sres)/sqrt(length(sres));

                    fprintf('iteration %d with residual %.2e and conditioning number %.2e\n',iter,nor/max([1,nor0])/obj.setting.newton.tol,con);
                end

                if retry==true
                    obj.unk.dotu.dof=dotu0;
                    obj.unk.u.dof=u0;
                    obj.unk.r.dof=r0;
                    obj.unk.rn.dof=rn0;
                    obj.unk.m.dof=m0;
                    obj.unk.dotpi.dof=dotpi0;
                    obj.unk.pi.dof=pi0;
                    obj.unk.al.dof=al0;
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

                obj.plot();

                fprintf("step completed\n");

                if good>=obj.setting.step.maxgood && obj.setting.step.lev+1<=obj.setting.step.nmax
                    obj.setting.step.lev=obj.setting.step.lev+1;
                end
            end

            fprintf("\nanalysis completed\n");

            close(obj.video);
        end
    end
end