classdef mylayerpot < handle
    %layerpotential claml.ss
    properties
        ss % aml.ssociated segments
        Dbdryp % double layer potential matrix at the boundary
        Dsbdryp % adjoint double layer potential matrix at the boundary
        modDsbdryp % adjoint double layer potential matrix at the boundary
        D % Kernel of D
        Ds %adjoint kernel
        S % single layer potential matrix at the boundary
        modS % modified single layer at the boundary
        dG % Gx+iGy evaluated at segdist
        G
        segdist
        density
        Mmat % M operator only for boundary
        M % # of segments
        tau
        type %'s' for single 'd' for double , 'ms' modified single, 'md' modified dstar
    end

    methods
        function ml = mylayerpot(ss, tau, type, phis,phitype)
            if nargin < 1, error('Segments required to define layer potential.'); end
            ml.ss = ss; ml.M = ss.M;
            ml.set_segdist;

            if nargin >= 2, ml.tau = tau; end
            if nargin >= 3, ml.type = type; end
            if nargin == 4
                phitype = [];
            end
            if nargin >= 4
               ml.assign_density(phis, phitype);
            end
        end

        
        % bdrys
        function dlayer_bdry(ml)
            if isempty(ml.Dbdryp) ==1 
                ml.Doperator;
                ml.Dbdryp = ml.D-0.5*eye(ml.ss.tN);  %interior boundary D- 0.5I;
            end
        end

        function dslayer_bdry(ml)
            if isempty(ml.Dsbdryp) ==1 
                ml.Dsoperator;
                ml.Dsbdryp = ml.Ds+0.5*eye(ml.ss.tN);%interior boundary D* + 0.5I;
            end
        end

        function modDslayer_bdry(ml)
            if isempty(ml.modDsbdryp) ==1 
                ml.Dsoperator;
                ml.setup_Mmat;
                IM = eye(size(ml.Ds))-ml.Mmat;
                ml.modDsbdryp = ml.Ds*IM+0.5*IM; %interior boundary D* + 0.5I modified;
                ml.Dsbdryp = ml.Ds+0.5*eye(ml.ss.tN);
            end
        end
        
        function slayer_bdry(ml)
                ml.Soperator;
        end

         function modSlayer_bdry(ml)
                ml.Soperator;
                ml.setup_Mmat;
                ml.modS =ml.S*(eye(ml.ss.tN)-ml.Mmat)+ml.Mmat; %modified single layer potential
         end

        % single and double layer evaluation not near bdry
        function [Dmat,deval] = dlayer(ml,tz)
            % Setup D[phi](z) matrix. multiplying by phi will give you
            % double layer potential
            if isempty(ml.density)
                error('density must be assigned to evaluate mylayerpot at bdry')
            end
            [tzr,tzc] = size(tz);
            tvN = tzr*tzc;
            tv = reshape(tz,[tvN,1]);
            dist = repmat(tv, [1 ml.ss.tN]) - repmat(ml.ss.zs.', [tvN 1]);
            [Gx,Gy] = myutils.gradG(dist,ml.tau);
            nuss = repmat(ml.ss.nus.', [tvN 1]);
            dnG = real(conj(nuss).*(-Gx-1i*Gy));
            Hs = repmat(ml.ss.ws,[tvN,1]);
            Dmat = dnG .* Hs; 
            if nargout <=1 
                return 
            end
            devalv = Dmat*ml.density;
            deval = reshape(devalv,[tzr,tzc]);
        end

        function [Smat,seval] = slayer(ml,tz)
            % Setup S[phi](z) matrix. multiplying by phi will give you
            % double layer potential
            if isempty(ml.density)
                error('density must be assigned to evaluate mylayerpot at bdry')
            end
            [tzr,tzc] = size(tz);
            tvN = tzr*tzc;
            tv = reshape(tz,[tvN,1]);
            dist = repmat(tv, [1 ml.ss.tN]) - repmat(ml.ss.zs.', [tvN 1]);
            Hs = repmat(ml.ss.ws,[tvN,1]);
            Smat = myutils.G(dist,ml.tau).* Hs; 
            if nargout ==1 
                return
            end
           sevalv = Smat*ml.density;
           seval = reshape(sevalv,[tzr,tzc]);
        end
        

         % setting up S, D, Ds
         function Doperator(ml)
            if isempty(ml.D) ==1
                ml.set_dG;
                nuss = repmat(ml.ss.nus.', [ml.ss.tN 1]);
                dnG = real(conj(nuss).*(-ml.dG));
                diagdnG = -ml.ss.kappas./(4*pi);
                dnG(myutils.dind(dnG)) = diagdnG;
                Hs = repmat(ml.ss.ws,[ml.ss.tN,1]);
                ml.D = dnG .* Hs; 
            end
        end

         function Dsoperator(ml)
            if isempty(ml.Ds) ==1
                ml.set_dG;
                nuss_star = repmat(ml.ss.nus, [1 ml.ss.tN]);
                dnGstar = real(conj(nuss_star).*(ml.dG));
                diagdnG = -ml.ss.kappas./(4*pi);
                dnGstar(myutils.dind(dnGstar)) = diagdnG; % adjoint have same diagonal.
                Hs = repmat(ml.ss.ws,[ml.ss.tN,1]);
                ml.Ds = dnGstar.*Hs;
            end
         end

         

         function Soperator(ml)
             if isempty(ml.S) ==1
                ml.set_G;
                sp = (ml.ss.speeds)./(2*pi);
                %prealloc
                circS2= zeros(ml.ss.tN); S1 = zeros(ml.ss.tN);
                %splitting near diagnoal
                for k = 1:ml.M
                    if ml.M ==1
                        currentN = ml.ss.tN;
                    else
                        currentN = ml.ss.Ns(k);
                    end
                    logmat = toeplitz(log(4*sin(pi*(0:currentN-1)/currentN).^2)); % create circulant matrix
                    tempS2 = (1/4/pi)*logmat;
                    tempS1 = toeplitz(quadr.kress_Rjn(currentN/2)).*(-1/4/pi);  % create circulant matrix
                    if ml.M==1
                        circS2 = tempS2;
                        S1 = tempS1;
                    else
                        circS2(ml.ss.indxs{k},ml.ss.indxs{k}) = tempS2;
                        S1(ml.ss.indxs{k},ml.ss.indxs{k}) = tempS1;
                    end
                end
                S2 = ml.G + circS2;
                v11 = myutils.dS2(ml.tau);
                S2(myutils.dind(S2)) = -log(sp*abs(v11))./2/pi;
                %consider weight for nystrom 
                S2 = S2.*repmat(ml.ss.ws,[ml.ss.tN 1]);
                S1 = S1.*repmat(sp, [ml.ss.tN 1]);
                ml.S = S1+S2;
                %testing S
                %Stild =  myutils.G(ml.segdist,ml.tau).*repmat(ml.ss.ws, [ml.ss.tN 1]);
                % compare S and Stild
                %errS = abs(ml.S-Stild); errS(dind(errS)) = 0; 
                %inferr = max(max(errS));
                %ml.S = (S1+S2).*repmat(ml.ss.ws, [ml.ss.tN 1]);
                %ml.S = S; 
             end
         end
        
        %green's function evaluations
        function set_dG(ml)
            if isempty(ml.dG) ==1
                ml.set_segdist;
                [Gx,Gy] = myutils.gradG(ml.segdist,ml.tau);
                ml.dG = Gx+1i*Gy;
            end
        end

        function set_G(ml)
            if isempty(ml.G) ==1
                ml.set_segdist;
                Geval = myutils.G(ml.segdist,ml.tau);
                ml.G =Geval;
            end
        end

        
         
        %helpers
        function setup_Mmat(ml)
            if isempty(ml.Mmat) == 1
                tN = ml.ss.tN;
                Med = zeros(tN);
                for i = 1:ml.M
                    s = ml.ss.segs{i};
                    Med(:,ml.ss.indxs{i}) = repmat(s.w,[tN 1]);
                end
                arclength = sum(ml.ss.ws);
                ml.Mmat = Med./arclength;
            end
        end
        
         function set_segdist(ml)
            if isempty(ml.segdist)==1  %reuse distance whenever possible 
                 dist = repmat(ml.ss.zs, [1 ml.ss.tN]) - repmat(ml.ss.zs.', [ml.ss.tN 1]);
                 ml.segdist = dist; 
           end
         end

        % evaluations
        function val = eval(ml, tz)
            if isempty(ml.density)
                error('density must be assigned to evaluate mylayerpot at bdry')
            end
            if nargin <= 2
            switch ml.type
                case {'s', 'ms'}
                    [~, val] = ml.slayer(tz);
                case 'd'
                    [~, val] = ml.dlayer(tz);
                otherwise
                    error('Invalid type: must be ''s'', ''d'', or ''ms''.');
            end
            end
        end

        function val = evalbd(ml,type)
            if isempty(ml.density)
                error('density must be assigned to evaluate mylayerpot at bdry')
            end
            if nargin >=2
                ml.type = type; % change type of mylayerpot;
            end
            switch ml.type
                case 's'
                    ml.slayer_bdry;
                    val = ml.S*ml.density;
                case 'd'
                    ml.dlayer_bdry;
                    val = ml.D*ml.density;
                case 'ms'
                    ml.modSlayer_bdry;
                    val = ml.modS*ml.density;
                otherwise
                    error('Invalid type: must be ''s'' or ''d''.');
            end
        end

      

        
        
        function density = process_density(ml, phis, phitype)
            % eval density depends on type
            if isempty(phitype)
                density = phis; 
                return
            end
            if ml.M~= numel(phis)
                error('there must be M many functions to evaluate for phi')
            end
            density = zeros(ml.ss.tN, 1);
            for i = 1:ml.M
                if phitype == 't'
                    density(ml.ss.indxs{i}) = phis{i}(2 * pi * ml.ss.ts(ml.ss.indxs{i}));
                elseif phitype == 'z'
                    density(ml.ss.indxs{i}) = phis{i}(ml.ss.zs(ml.ss.indxs{i}));
                end
            end
        end
        

        %related to class setup and update
        function assign_density(ml,phis,phitype)
            if isempty(ml.density)
                if nargin <=2
                    phitype = [];
                end
                ml.density = ml.process_density(phis,phitype);
            end
        end
        
        function requadrature(ml,Ns)
            ml.reset; % reset everything that depends on tau
            ml.ss.requadrature(Ns);
            ml.segdist = repmat(ml.ss.zs, [1 ml.ss.tN]) - repmat(ml.ss.zs.', [ml.ss.tN 1]);
            ml.Mmat=[]; % M operator 
            ml.density = []; % size of density also changes
        end

        function update_density(ml,phis,phitype)
            ml.density = [];
            if nargin <=2
                phitype = [];
            end
            ml.assign_density(phis,phitype);
        end

        function updatetau(ml,tau)
            ml.reset;
            ml.tau = tau;
        end

        function reset(ml) %reset everything else that depends on tau
            ml.Dbdryp = []; 
            ml.Dsbdryp = []; 
            ml.modDsbdryp = [];
            ml.modS = [];
            ml.D = [];
            ml.Ds = []; %adjoint kernel
            ml.S = []; % single layer potential matrix at the boundary
            ml.G = [];
            ml.dG = [];
        end
    end
end