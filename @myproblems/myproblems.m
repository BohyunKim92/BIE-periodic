classdef myproblems < handle
    %PROBLEMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ss % segmets used
        ml % layerpotential class
        type % 'D' for dirichlet, 'N' for neumann , 'S' for steklov
        A % double layer at bdry
        B % used for steklov evp
        rhs % only used for dirchlet and neumann
        co %corresponds to density + fluxes only used for dirichlet and neumann 
        phis % recovered density
        evals %only used for steklov evp stores eigenvalues
        efuncs %only used for steklov evp stores eigenfunctions
        fluxes % calculated flux
        %errsinf % error for current segments with infnorms
        %abserr % error for current segment with abserr for each component
        %errres % residual error for solving the density
        %converrs %multiple errors varying size of N
        gtype % type of bdry. t takes s.t, z takes s.z...etc
        gs %the function for boundary
        M %# of segments
        Ns
        tau
    end
    
    methods
        function pr = myproblems(mlp,type)
            if nargin==0, return; end % empty constructor (for copy)
            pr.ss = mlp.ss;
            pr.ml = mlp;
            pr.M = mlp.M;
            pr.Ns = mlp.ss.Ns;
            pr.tau = mlp.tau;
            if nargin ==2
                pr.type = type; return;
            end
        end

        function solve(pr)
            setup_A(pr);
            if pr.type == 'S'
                [pr.efuncs,D] = eig(pr.A,pr.B);
                pr.evals = diag(D);
            elseif pr.type =='D'
                % dirichlet bvp
                if isempty(pr.rhs) || isempty(pr.A)
                    error('make sure to set up the rhs and matrix solver first.')
                end
                pr.co = pr.A\pr.rhs;
                pr.phis = pr.co(1:pr.ss.tN);
                pr.ml.update_density(pr.phis);
                pr.fluxes = 0;
                if pr.M ~=1
                    pr.fluxes = pr.co(pr.ss.tN+1:pr.ss.tN+pr.M);
                end

            elseif pr.type =='N'
                %neumann bvp
                if isempty(pr.rhs) 
                    error('make sure to set up the rhs and matrix solver first.')
                end
                pr.co = pr.A\pr.rhs;
                pr.phis = pr.co;
                pr.fluxes = zeros(pr.M,1);
                for i = 1:pr.M
                    pr.fluxes(i) = pr.ss.ws(pr.ss.indxs{i})* pr.phis(pr.ss.indxs{i});
                end
                
                %steklov evp
            else
                error(['the problem must be specified by D,N,S']);
            end
        end
        
        function setuprhs(pr,fs,gtype)
            if nargin <=2 
                if isa(class(fs),'mylayerpot')
                    val = fs.eval;
                else
                    val = fs;
                end
                if pr.type =='D' && pr.M >=2
                    pr.rhs = zeros(pr.ss.tN+pr.M,1);
                    pr.rhs(1:pr.ss.tN) = val;
                else
                    pr.rhs = val;
                end
                return
            end
            if  numel(fs) ~= pr.M
                error(['there must be ', num2str(pr.M),' functions for the boundary condition'])
            end
            pr.gs = fs;
            pr.gtype = gtype;
            if ~(pr.type == 'D'|| pr.type ~= 'N'|| pr.type ~= 'S') 
                error('you must setup correct problem type');
            end
            if gtype == 't'
                 target = 2*pi*pr.ss.ts;
            elseif gtype == 'z'
                target = pr.ss.zs;
            end
            pr.rhs = myutils.evalf(fs,target,pr.ss.indxs);
        end
        
       
       function setup_A(pr)
           if isempty(pr.M)
               error('must provide at least one segment to build bdry equation using double layer potential!');
           end
           if pr.type == 'D'
               if pr.M ==1
                   pr.ml.dlayer_bdry; 
                   pr.A = pr.ml.Dbdryp; %assign
               else
                   pr.dirichlet_mult; % use double layer bdry equations with \sum AjG(z-\beta_j) form
               end
           elseif pr.type == 'N'
               pr.ml.dslayer_bdry
               pr.A = pr.ml.Dsbdryp;
           elseif pr.type == 'S'
               pr.ml.modDslayer_bdry;
               pr.ml.modSlayer_bdry;
               pr.B = pr.ml.modS;
               pr.A = pr.ml.modDsbdryp;
           else
               error('Only Dirichlet, Neumann, and Steklov Supported');
           end
       end
      
         

         function CDbdrys = cdlayer_bdry_single(pr)
            %set up complex layer potential
            % dlogtheta = myutils.dlogtheta(dist,ss.tau);
            % zpss = repmat((ss.zps./(2*pi)).', [ss.tN 1]);
            % dnG = real((1/(2*pi*1i))*(zpss.*dlogtheta));
            % sp = ss.speeds./(2*pi);
            % dia = -(ss.kappas./(4*pi)).*sp;
            % dnG(myutils.dind(dnG)) = dia;
            CDbdrys = [];
         end
        

        function dirichlet_mult(pr)
        % assume multiple holes set up bdry eqn
        % fill all block matrices where Dbrym = [D+0.5I, B; C, E]
      
        Dbdrym  = zeros(pr.ss.tN+pr.M);
        %fill C
        C = zeros(pr.ss.M,pr.ss.tN);
        for i = 1:pr.ss.M
            if i <=pr.ss.M-1
                C(i, pr.ss.indxs{i}) =pr.ss.ws(pr.ss.indxs{i});
            end
        end
        
        %fill E
        E = zeros(pr.ss.M); E(pr.ss.M,:) = 1;
        
        %fill B
        dist1 = repmat(pr.ss.zs,[1,pr.ss.M])-repmat((pr.ss.cs).',[pr.ss.tN,1]);
        B = myutils.G(dist1,pr.tau);
        
        %fill D
        % if isempty(pr.ml.Dbdryp) ==1
        %     pr.ml.dlayer_bdry; % calculate  only when empty
        % end
        pr.ml.dlayer_bdry;
        Dbdrym(1:pr.ss.tN,1:pr.ss.tN) = pr.ml.Dbdryp;
        Dbdrym(1:pr.ss.tN,pr.ss.tN+1:pr.ss.tN+pr.ss.M) = B;
        Dbdrym(pr.ss.tN+1:pr.ss.tN+pr.ss.M,1:pr.ss.tN) = C; 
        Dbdrym(pr.ss.tN+1:pr.ss.tN+pr.ss.M,pr.ss.tN+1:pr.ss.tN+pr.ss.M) = E;
        pr.A = Dbdrym;
        end
    
        function uu = eval(pr,zz,eigsets)
            if nargin >=3
                if pr.type ~='S'
                    error('must be steklov evp')
                end
            end
            if pr.type =='D'
                [~, uu] = pr.ml.dlayer(zz);
                if pr.M== 1
                    return
                else
                    %must consider fluxes for M>1 case
                     [tzr,tzc] = size(zz);
                     tvN = tzr*tzc;
                     tv = reshape(zz,[tvN,1]);
                     % u1vec = Dmat*pr.phis;
                     dist = repmat(tv,[1,pr.M])- repmat(pr.ml.ss.cs.',[tvN,1]);
                     u2vec =  myutils.G(dist,pr.tau)*pr.fluxes;
                     uu2 = reshape(u2vec,[tzr,tzc]);
                    uu = uu+uu2;
                end
            elseif pr.type == 'N'
                % the final solution is 
                [~, uu] = pr.ml.slayer(zz,pr.phis);
            end
        end
        
        function eval = calc_stream1(pr,z0,tend) % increasing in z
            sol = ode23(@(t,z) pr.dzdt(t,z),[0 tend],z0);  
            eval= deval(sol,linspace(0,max(sol.x),50));
        end

        function eval = calc_stream2(pr,z0,tend) % decreasing in z
            sol = ode23(@(t,z) -pr.dzdt(t,z),[0 tend],z0);  
            eval= deval(sol,linspace(0,max(sol.x),50));
        end

        function plot_stream1(pr,target,tend,ptau)   
        for i = 1:numel(target)
            z0 = target(i);
            eval = calc_stream1(pr,z0,tend);
            ind = ptau.intau(eval);
            plot(eval(ind),'k','LineWidth',2); hold on;
            %hold on; plot(eval(1),'k*');
            % notind = ~ind; 
            % noteval = eval(notind);
            % if sum(notind)>0 % if there is a segment that is outside of the Ptau
            %     shifted = ptau.all_dir_shift(noteval);
            %     for k = 1:numel(shifted)
            %         if isempty(shifted{k}) == 0
            %             plot(shifted{k},'r','LineWidth',2); hold on;
            %         end
            %     end
            % end
        end
        end

        function plot_stream2(pr,target,tend,ptau)   
        for i = 1:numel(target)
            z0 = target(i);
            eval = calc_stream2(pr,z0,tend);
            ind = ptau.intau(eval);
            plot(eval(ind),'k','LineWidth',2); hold on;
        end
        end

        function g = dzdt(pr,t,z)
            u = pr.eval(z);
            xu = pr.eval(z+0.001);
            yu = pr.eval(z+0.001*1i);
            gradux = (xu-u)/0.001;
            graduy = (yu-u)/0.001;
            absg = abs(gradux+1i*graduy);
            g = (gradux + 1i*graduy)./absg;
        end
       

        function plot(pr,pt,levels,color)
            % plot solution given the target. 
            if pr.type == 'D'||pr.type == 'N'
                if isempty(pt.solgrid)
                    pt.add_grid;
                end
                zz = pt.solgrid;
                uu = pr.eval(zz);
                ind = intau(pt,pt.solgrid,2);
                uu(ind) = NaN;
                lw = 2; %line width of the contour
                figure;
                xx = real(zz);
                yy = imag(zz);
                lw = 2;
                contour(xx,yy,uu,levels,'LineWidth',lw,'color',color);hold on;
                % fill in hole
                pr.ss.plot2(2); hold on;
                % plot torus boundary
                pt.plot2(1);
                axis equal;
                ax = gca;
                set(gca,'xtick',[],'ytick',[]);
                set(gca,'fontsize',20)
            end
        end
         
        function requadrature(pr, Ns)
            %change the number of grid being used
            pr.Ns = Ns;
            pr.ml.requadrature(Ns); %wipe everything
            pr.A =[]; % double layer at bdry
            pr.B  = []; % used for steklov evp
            pr.rhs = []; % only used for dirchlet and neumann
            pr.co=[]; %corresponds to density + fluxes only used for dirichlet and neumann 
            pr.phis=[]; % recovered density
        pr.evals=[]; %only used for steklov evp stores eigenvalues
        pr.efuncs=[]; %only used for steklov evp stores eigenfunctions
        pr.fluxes= []; 
        end
        
        function updatetau(pr,tau)
            pr.tau = tau;
            pr.ml.updatetau(tau);
            pr.A =[]; % double layer at bdry
            pr.B  = []; % used for steklov evp
            pr.rhs = []; % only used for dirchlet and neumann
            pr.co=[]; %corresponds to density + fluxes only used for dirichlet and neumann 
            pr.phis=[]; % recovered density
            pr.evals=[];%only used for steklov evp stores eigenvalues
            pr.efuncs=[]; %only used for steklov evp stores eigenfunctions
            pr.fluxes= []; 
        end
    end
    
end