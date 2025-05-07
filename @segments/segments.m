classdef segments < handle
    %SEGMENTS Class for handling multiple boundary segments in 2D.
    %   This class manages geometry, quadrature information, and plotting 
    %   utilities for a collection of parametric boundary segments.
    
    properties
        segs        % Cell array of segment structures
        orient      % 'c' for clockwise, 'cc' for counterclockwise
        Ns          % Number of quadrature nodes per segment
        zs          % All quadrature points combined (complex)
        ws          % Quadrature weights (row vector)
        ts          % Parametrization values in [0,1] (column vector)
        speeds      % |dZ/dt| at quadrature points
        kappas      % Curvature at quadrature points (optional)
        nus         % Outward normal vectors at quadrature points
        M           % Number of segments
        indxs       % Cell array of index sets for each segment
        tN          % Total number of quadrature points
        cs          % Center points of segments
        rs          % Radii of segments
        zps         % Derivatives dZ/dt at quadrature points
        arclength   % Total arclength (sum of quadrature weights)
        pltzs       % Quadrature points for plotting
        pltnus      % Normals for plotting
        pltind      % Index sets for plotting
        oscils      % Wavenumbers for star-shaped oscillating boundaries
    end
    
    methods
        function ss = segments(segms)
            % Constructor: Initialize segment object from a cell array
            if nargin == 0, return; end
            M = size(segms, 1); ss.M = M;
            Ns = zeros(M,1); segs = cell(M,1);

            % Convert counter-clockwise to clockwise if needed
            if (isempty(ss.orient) || ss.orient == 'cc')
                for i = 1:M
                    s = segms{i};
                    s.x = flipud(s.x);
                    s.t = flipud(s.t);
                    s.nx = -flipud(s.nx);
                    s.w = flipud(s.w);
                    s.kappa = -flipud(s.kappa);
                    s.speed = flipud(s.speed);
                    Ns(i) = size(s.x,1);
                    segs{i} = s;
                end
                ss.orient = 'c'; % now oriented clockwise
            end

            ss.Ns = Ns;
            ss.tN = sum(Ns);
            ss.segs = segs;
            ss.indxs = ss.indices;
            ss.set_vectors;
        end

        function set_vectors(ss)
            % Assigns all vector properties from individual segments
            ss.zs = zeros(ss.tN,1); ss.ws = zeros(1,ss.tN);
            ss.ts = zeros(ss.tN,1); ss.speeds = zeros(1,ss.tN);
            ss.nus = zeros(ss.tN,1); ss.zps = zeros(ss.tN,1);
            ss.kappas = zeros(ss.tN,1);

            for i = 1:ss.M
                s = ss.segs{i};
                if ~isempty(s.Zp)
                    ss.zps(ss.indxs{i}) = -s.Zp(s.t);
                end
                ss.zs(ss.indxs{i}) = s.x;
                ss.nus(ss.indxs{i}) = s.nx;
                ss.ws(ss.indxs{i}) = s.w;
                ss.ts(ss.indxs{i}) = s.t;
                ss.speeds(ss.indxs{i}) = s.speed;
                ss.kappas(ss.indxs{i}) = s.kappa;
            end

            ss.arclength = sum(ss.ws);
            if ~isempty(ss.pltzs) || ~isempty(ss.pltnus) || ~isempty(ss.pltind)
                ss.pltzs = [];
                ss.pltnus = [];
                ss.pltind = [];
                ss.setpltvars;
            end
        end

        function requadrature(ss, nNs)
            % Change quadrature resolution for each segment
            if size(nNs) ~= [ss.M, 1]
                error('nNs must be a column vector of size [M,1].');
            end
            ss.Ns = nNs; ss.tN = sum(nNs);
            ss.indxs = ss.indices;

            for i = 1:ss.M
                s = ss.segs{i};
                s.requadrature(nNs(i));
                s.x = flipud(s.x);
                s.t = flipud(s.t);
                s.nx = -flipud(s.nx);
                s.w = flipud(s.w);
                s.kappa = -flipud(s.kappa);
                s.speed = flipud(s.speed);
                ss.segs{i} = s;
            end

            ss.set_vectors;
            ss.arclength = sum(ss.ws);
        end

        function newss = copy(ss)
            % Create a deep copy of the segments object
            newss = segments();
            newss.M = ss.M;
            newss.segs = cellfun(@myutils.copy_segment, ss.segs, 'UniformOutput', false);
            newss.indxs = ss.indxs;
            props = {'orient','Ns','zs','ws','ts','speeds','kappas','nus',...
                     'tN','cs','rs','zps','arclength','pltzs','pltnus','pltind','oscils'};
            for k = 1:length(props)
                newss.(props{k}) = ss.(props{k});
            end
        end

        function indxs = indices(ss, Ns)
            % Return index sets for each segment
            if nargin < 2, Ns = ss.Ns; end
            indxs = cell(length(Ns),1);
            for i = 1:length(Ns)
                if i == 1
                    indxs{i} = 1:Ns(i);
                else
                    indxs{i} = (1+sum(Ns(1:i-1))):sum(Ns(1:i));
                end
            end
        end

        function [zzind] = inseg(ss, zz, opt)
            % Determine if points zz are inside (opt=1) or outside (opt=2) the segments
            if nargin < 3, opt = 1; end
            ind = false(size(zz));
            for i = 1:ss.M
                s = ss.segs{i};
                tval = myutils.ztot(zz, ss.cs(i));
                rseg = abs(s.Z(tval) - ss.cs(i));
                rz = abs(zz - ss.cs(i));
                ind(rz <= rseg) = true;
            end
            zzind = opt == 1 && ind | opt == 2 && ~ind;
        end

        function plot(ss)
            % Plot segments with tangents and normal vectors
            figure;
            for i = 1:ss.M
                s = ss.segs{i};
                cx = real(s.x); cy = imag(s.x);
                cdiff = circshift(s.x,-1) - s.x;
                quiver(cx, cy, real(cdiff), imag(cdiff), 0, 'k-', 'LineWidth', 2); hold on;
                quiver(cx(1), cy(1), real(cdiff(1)), imag(cdiff(1)), 0, 'r', 'LineWidth', 2);
                quiver(cx, cy, real(s.nx), imag(s.nx), 'b-');
                if ~isempty(ss.zps)
                    quiver(cx, cy, real(ss.zps(ss.indxs{i})), imag(ss.zps(ss.indxs{i})), 'r-', 'LineWidth', 2);
                end
            end
            h = arrayfun(@(i) plot(real(ss.cs(i)), imag(ss.cs(i)), '*', 'DisplayName', num2str(ss.cs(i))), ...
                         1:ss.M, 'UniformOutput', false);
            legend([h{:}], 'Location', 'best');
            title('Segments with normals and tangents');
            hold off;
        end

        function plot2(ss, opt, width, adjust1, adjust2)
            % Pretty plotting for publication, with gaps or fills
            if nargin < 2, opt = 1; end
            if nargin < 3, width = 2.5; end
            if nargin < 4, adjust1 = 0.01; end
            if nargin < 5, adjust2 = 0.017; end

            ss.setpltvars;
            segs = ss.pltzs - adjust1 * ss.pltnus;
            wsegs = ss.pltzs - adjust2 * ss.pltnus;

            for i = 1:ss.M
                ith = ss.pltind{i};
                if opt == 1
                    plot(wsegs(ith), 'w', 'LineWidth', width); hold on;
                    plot(segs(ith), 'color',[.7 .7 1], 'LineWidth', width); hold on;
                elseif opt == 2
                    fill(real(wsegs(ith)), imag(wsegs(ith)), 'w'); hold on;
                    plot(wsegs(ith), 'w', 'LineWidth', width); hold on;
                    fill(real(segs(ith)), imag(segs(ith)), [.7 .7 1]); hold on;
                    plot(segs(ith),'color', [.7 .7 1], 'LineWidth', width); hold on;
                end
            end
        end

        function setpltvars(ss)
            % Precompute variables for plotting with loops
            segs = zeros(ss.tN + ss.M, 1);
            nus = zeros(ss.tN + ss.M, 1);
            Ns_plot = ss.Ns + 1;
            indxp = ss.indices(Ns_plot);
            for i = 1:ss.M
                idx = ss.indxs{i};
                segs(indxp{i}(1:ss.Ns(i))) = ss.zs(idx);
                segs(indxp{i}(end)) = ss.zs(idx(1));
                nus(indxp{i}(1:ss.Ns(i))) = ss.nus(idx);
                nus(indxp{i}(end)) = ss.nus(idx(1));
            end
            ss.pltzs = segs;
            ss.pltnus = nus;
            ss.pltind = indxp;
        end

        function centers(ss, cs)
            ss.cs = cs;
        end

        function radiuss(ss, rs)
            ss.rs = rs;
        end

        function oscil(ss, oscil)
            ss.oscils = oscil;
        end
    end

    methods(Static)
        [ss] = circles(cs, rs)
        [ss] = trefs(cs, rs)
        [ss] = neumann_regions
        [ss] = steklov_many_holes
        s = mycircle(c, r)
        s = mystar(c, r, osc)
        s = mytref(c, r)
        indxs = coarse_to_fine_indxs(finess, coarsess)
    end
end
