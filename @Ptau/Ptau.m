classdef Ptau < handle
    % Ptau stores information about the torus and fundamental domain

    properties
        tau         % The period information (complex number)
        slp         % Slope: imag(tau)/real(tau)
        fs          % Function used to determine the interior of the fundamental domain
        shift_fs    % Shifting functions: rd (right down), ru (right up), etc.
        shift_dir   % Shifting directions: rd (right down = +1 - tau), ru (right up), etc.
        offsetx     % X-offset for coordinate shifting
        indxs       % Indices of each segment with respect to the total number of quadrature points
        tN          % Total number of quadrature points for generating the fundamental domain
        fsbdry      % Cell array of functions: {lf, rf, bf, tf}, which output the boundary curves
        zs          % Coordinates combined into z, used for the boundary
        ts          % Quadrature parameter values in [0, 1] (column vector) used for boundary
        solgrid     % Coordinates for the solution inside the fundamental domain.
        Ns          % Number of quadrature points for each line segment of the parallelogram
        M           % Number of segments (in this case, 5)
    end
    
    methods
    function pt = Ptau(tau, Ns)
        % Construct an instance of this class
        if nargin == 0
            return;  % Empty constructor (for copying)
        end
        pt.M = 4;
        if nargin == 1
            Ns = 100 * ones(pt.M, 1); % default number of quadrature
        end
        pt.tau = tau;
        pt.slp = imag(tau) / real(tau);
        pt.offsetx = -real(pt.tau) / 2;
        pt.Ns = Ns;
        pt.tN = sum(pt.Ns);
        pt.indxs = pt.indices;
        pt.set_bdry_functions;
        pt.set_vectors;
    end

    function set_bdry_functions(pt)
        % Assign boundary functions for fs (for plotting purposes)
        if isempty(pt.fsbdry)  % Ensure calculation occurs only once
            if real(pt.tau) == 0
                lf = @(y) imag(pt.tau) * 1i * y;
                rf = @(y) 1 + imag(pt.tau) * 1i * y;
                bf = @(x) x;
                tf = @(x) x + pt.tau;
            else
                yintc2 = pt.slp * (-1);
                lf = @(x) (x + 1i * (pt.slp * x) + pt.offsetx);
                rf = @(x) (x + 1i * (pt.slp * x + yintc2) + pt.offsetx);
                bf = @(x) x + pt.offsetx;
                tf = @(x) x + pt.offsetx + 1i * imag(pt.tau);
            end
            pt.fsbdry = {bf, rf, tf, lf};
        end
    end

    function set_intau_functions(pt)
        if isempty(pt.fs) || isempty(pt.shift_fs) || isempty(pt.shift_dir)
            if real(pt.tau) == 0
                lf = @(z) real(z) > 0;
                rf = @(z) real(z) < 1;
                bf = @(z) imag(z) > 0;
                tf = @(z) imag(z) < imag(pt.tau);
            else
                yintc2 = pt.slp * (-1);
                lf = @(z) imag(z) < pt.slp * real(z - pt.offsetx);
                rf = @(z) imag(z) > pt.slp * real(z - pt.offsetx) + yintc2;
                bf = @(z) imag(z) > 0;
                tf = @(z) imag(z) < imag(pt.tau);
            end
            pt.fs = {bf, rf, tf, lf};

            % Create shifting function
            pt.shift_fs = {
                @(z) ~lf(z) & ~tf(z);  % Left up position -> require rd shifting
                @(z) ~lf(z) & ~bf(z);  % Left bottom -> require ru shifting
                @(z) ~lf(z) & bf(z) & tf(z);  % Right side
                @(z) ~rf(z) & ~bf(z);  % Left up side
                @(z) ~rf(z) & bf(z) & tf(z);  % Left side
                @(z) ~rf(z) & ~tf(z);  % Left down side
                @(z) rf(z) & lf(z) & ~bf(z);  % Top side
                @(z) rf(z) & lf(z) & ~tf(z);  % Bottom side
            };

            % Corresponding shifting directions
            pt.shift_dir = {
                1 - pt.tau;  % rd
                1 + pt.tau;  % ru
                1;           % Right side
                -1 + pt.tau; % lu
                -1;           % Left side
                -1 - pt.tau;  % ld
                pt.tau;       % Top side
                -pt.tau       % Bottom side
            };
        end
    end

    function [zzind] = intau(pt, zz, opt)
        % Given coordinates zz, if opt = 1, return index of zz in P_tau
        % If opt = 2, return index of zz outside P_tau (P_tau = parallelogram)
        pt.set_intau_functions;
        if nargin == 1
            error('intau requires zz');
        end
        if nargin <= 2
            opt = 1;
        end

        zzind = logical(pt.fs{1}(zz) .* pt.fs{2}(zz) .* pt.fs{3}(zz) .* pt.fs{4}(zz));
        if opt == 2
            zzind = ~zzind;  % Get the outside of the Ptau
        elseif opt ~= 2 && opt ~= 1
            error('opt must be 1 or 2');
        end
    end

    function set_vectors(pt)
        if isempty(pt.zs) || isempty(pt.ts)  % Only calculate once
            % Initialize
            pt.zs = zeros(pt.tN, 1);
            pt.ts = zeros(pt.tN, 1);
            for i = 1:pt.M
                % Clockwise orientation
                if real(pt.tau) == 0
                    if i <= 2
                        ti = linspace(0, 1, pt.Ns(i));
                    else
                        ti = linspace(1, 0, pt.Ns(i));
                    end
                else
                    if i == 1
                        ti = linspace(0, 1, pt.Ns(i));
                    elseif i == 2
                        ti = linspace(1, 1 + real(pt.tau), pt.Ns(i));
                    elseif i == 3
                        ti = linspace(1 + real(pt.tau), real(pt.tau), pt.Ns(i));
                    elseif i == 4
                        ti = linspace(real(pt.tau), 0, pt.Ns(i));
                    end
                end
                pt.ts(pt.indxs{i}) = ti;
                pt.zs(pt.indxs{i}) = pt.fsbdry{i}(ti);
            end
        end
    end

    function add_grid(pt,xcount,ycount)
        % create grid for plotting solution 
        if nargin <=1
            xcount = 100; % default grid resolution to plot 
            ycount = 100;
        elseif nargin == 2
            ycount = xcount;
        end
        xlims = [pt.offsetx+0.01, 1+real(pt.tau)+pt.offsetx-0.01];
        ylims = [0, imag(pt.tau)];
        x = linspace(xlims(1),xlims(2),xcount); y = linspace(ylims(1),ylims(2),ycount); 
        [xx,yy] = meshgrid(x,y); pt.solgrid = xx+1i*yy;
    end

    function indxs = indices(pt)
        % Suppose Ns = [N1, N2, ..., NM]. Then indices for each segment.
        indxs = cell(pt.M, 1);
        for i = 1:pt.M
            if i == 1
                indxs{i} = 1:pt.Ns(i);
            else
                indxs{i} = 1 + sum(pt.Ns(1:i - 1)):sum(pt.Ns(1:i));
            end
        end
    end

    function shiftedzz = all_dir_shift(pt, zz)
        n = numel(pt.shift_fs);
        shiftedzz = cell(n, 1);
        pt.set_intau_functions;
        for k = 1:n
            shiftedzz{k} = pt.directional_shift(zz, k);
        end
    end

    function result = directional_shift(pt, zz, opt)
        % Given zz (the location), apply shifting with option opt
        % opt 1 = rd shift, opt 2 = ru shift, and so on (see shift_fs)
        pt.set_intau_functions;
        if opt < 1 || opt > numel(pt.shift_fs)
            error('Invalid option for shifting.');
        end
        % Get corresponding mask function and shift vector
        shiftFunc = pt.shift_fs{opt};
        shiftdir = pt.shift_dir{opt};
        ind = shiftFunc(zz);
        result = unique(zz(ind) + shiftdir);
    end

    function plot(pt, width)
    % Plot the boundary of the fundamental domain along with the tangent
    % vectors.
        if nargin <= 1
            width = 4; % Set the line thickness for the boundary
        end
        pt.set_bdry_functions;
        pt.set_vectors;
        for i = 1:pt.M
            cx = real(pt.zs(pt.indxs{i}));
            cy = imag(pt.zs(pt.indxs{i}));
            plot(cx, cy, 'k', 'LineWidth', width, 'MarkerSize', 2);
            hold on;
            mid = floor(numel(cx) / 2);
            fct = floor(pt.Ns(i) / 10);
            quiver(cx(mid - fct), cy(mid - fct), cx(mid + fct) - cx(mid - fct), cy(mid + fct) - cy(mid - fct), 0, 'r', 'LineWidth', width, 'MarkerSize', 2);
            hold on;
        end
        title('The Fundamental Domain');
    end

    function plot2(pt, opt, width)
        % Used to clean near the boundary of the segment to generate a pretty plot
        % opt == 1: plot the boundary only
        % opt == 2: plot the boundary and shade inside

        pt.set_bdry_functions;
        pt.set_vectors;
        % Parameter setup
        if nargin <= 1
            opt = 1;
            width = 4; % thickness of the line for boundary
        elseif nargin == 2
            width = 4;
        end

        if opt == 1
            for i = 1:pt.M
                cx = real(pt.zs(pt.indxs{i}));
                cy = imag(pt.zs(pt.indxs{i}));
                plot(cx, cy, 'r', 'LineWidth', width, 'MarkerSize', 2);
                hold on;
            end
        elseif opt == 2
            x = [pt.offsetx, 1 + pt.offsetx, 1 + real(pt.tau) + pt.offsetx, real(pt.tau) + pt.offsetx];
            y = [0, 0, imag(pt.tau), imag(pt.tau)];
            fill(x, y, [.5 .5 .5], 'EdgeColor', 'none');
            alpha(0.2);
            hold on;
            for i = 1:pt.M
                cx = real(pt.zs(pt.indxs{i}));
                cy = imag(pt.zs(pt.indxs{i}));
                plot(cx, cy, 'r', 'LineWidth', width, 'MarkerSize', 2);
                hold on;
            end
        end
        axis equal;
    end
    end
end
