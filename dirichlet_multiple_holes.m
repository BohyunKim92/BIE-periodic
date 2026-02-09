% Dirichlet Problem: Multiple Holes Example
% Solves Laplace's equation with Dirichlet boundary conditions
% Example 2 in Section 4 of 
% "Layer Potential Methods for Doubly-Periodic Harmonic Functions"
%
% This script uses the following classes/functions:
% - @segment, @quadr, @pointset, and selected utilities from @util in the 
%   mpspack library (by Alex Barnett and Timo Betcke, 2008–2012)
% - Custom methods written by Bohyun Kim (University of Utah, Math, 2025)

clear all; close all; clc;
addpath('mps_helpers');  % Add helper functions to path
%addpath([pwd,'/data']);

%% Define the geometry: three circular holes(casenum = 1), three trefoils (all other cases)
casenum = 2
M = 3;
Ns = [50; 50; 50]; % quadrature points
cs = [0.7+0.5*1i;0.3+ 0.3*1i;0]; % centers
rs = [0.1; 0.1; 0.1]; % radius'

% Define periodicity (tau = complex lattice vector)
%tau = 1i                     % Option: square torus
tau = 0.5 + sqrt(3)/2 * 1i;    % Option: equilateral torus

%define boundary
As = [3;-1;-2];
fexact = @(z) As(1)*myutils.G(z-cs(1),tau)+As(2)*myutils.G(z-cs(2),tau)+As(3)*myutils.G(z-cs(3),tau);
fs = {fexact;fexact;fexact};
rphi = {@(t) -10*sin(t); @(t) 10* sin(3*t); @(t) -10* sin(t)};

if casenum ==1
    ss = segments.circles(cs, rs);          % Create circular segment
    ss2 = segments.circles(cs,[0.18,0.18,0.18]);
else
    ss = segments.trefs(cs, rs);
    ss2 = segments.trefs(cs,[0.18,0.18,0.18]);
end
ss.requadrature(Ns)
ss2.requadrature([1250; 1250; 1250])

%% Build layer potential representation
ml = mylayerpot(ss, tau,'d');                         % define layer potentials
rhss = ss.copy;
rhs = mylayerpot(rhss, tau, 's', rphi, 't');   
rhseval = rhs.evalbd;   
rhseval = rhseval+fexact(rhss.zs);        % evaluate boundary    
% Solve Dirichlet problem
pr = myproblems(ml, 'D');                  % 'D' = Dirichlet
pr.setuprhs(rhseval);                      % Set RHS
pr.solve;                                  % Solve linear system

pt = Ptau(tau); % set up parallelogram domain
xcount = 100; ycount = 100;
pt.add_grid(xcount,ycount)
levels = As(1)*[-0.2:0.1:0.2,-0.01:0.05:0.1,0.16, -0.05:0.02:0.05];
pr.plot(pt,levels)


%% Error calculation
uexact = rhs.eval(ss2.zs)+fexact(ss2.zs);
u = pr.eval(ss2.zs);
error = norm(u-uexact,'inf')