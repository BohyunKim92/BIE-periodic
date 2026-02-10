% Neumann Problem: Multiple Hole Example
% Solves Laplace's equation with Neumann boundary conditions
% Example 3 in Section 4 of 
% "Layer Potential Methods for Doubly-Periodic Harmonic Functions"
%
% This script uses the following classes/functions:
% - @segment, @quadr, @pointset, and selected utilities from @util in the 
%   mpspack library (by Alex Barnett and Timo Betcke, 2008–2012)
% - Custom methods written by Bohyun Kim (University of Utah, Math, 2025)

clear all; close all; clc; verb = 0;
addpath('mps_helpers');  % Add helper functions to path

%% Define geometry: Neumann region
ss = segments.neumann_regions; %segment defined in @segments folder
Ns = 50*ones(ss.M,1); % Coarse discretization
ss.requadrature(Ns);
as = [1,2,3,4,-5,-2,-3,0];
% Define periodicity (tau = complex lattice vector)
%tau = 1i;                     % Option: square torus
tau = 0.5 + sqrt(3)/2 * 1i;    % Option: equilateral torus


%% Build layer potential representation
ml = mylayerpot(ss,tau,'d');

% Solve Neumann problem
pr = myproblems(ml,'N');
rhs = rhs1(ss,as,tau);
pr.setuprhs(rhs);
pr.solve;

pt = Ptau(tau); % set up parallelogram domain
xcount = 100; ycount = 100;
pt.add_grid(xcount,ycount)
levels = [linspace(-5,4,45),1e-6];
pr.plot(pt,levels)

%%
%compute error at the boundary
ml2= mylayerpot(ss,tau,'s',pr.phis); % define modified single layer to calculate
u_bdry = ml2.evalbd;
uexact = uexact1(ss.zs,ss,as,tau); % exact solution
err = u_bdry-uexact;
m = mean(err);
new = u_bdry-m;
err2 = new-uexact;
inferr = norm(err2,'inf');
fprintf('The calculated error at the boundary is %d\n', inferr);

function rhs = rhs1(ss,as,tau)
% test neumann bvp data
zshifts = repmat(ss.zs,[1,ss.M])- repmat(ss.cs.',[ss.tN, 1]);
[Gx,Gy] = myutils.gradG(zshifts,tau);
g = real(conj(repmat(ss.nus,[1,ss.M])).*(Gx+1i*Gy)); %dG/dnu(z)
rhs = g*as.';
end

function u = uexact1(tz,ss,as,tau)
[tzr,tzc] = size(tz);
tvN = tzr*tzc;
tv = reshape(tz,[tvN,1]);
zshifts= repmat(tv, [1,ss.M]) - repmat(ss.cs.', [tvN 1]);
Geval = myutils.G(zshifts,tau);
uvec = Geval*as.';
u = reshape(uvec,[tzr,tzc]); 
end
