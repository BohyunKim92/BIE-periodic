% Dirichlet Problem: Single Hole Example
% Solves Laplace's equation with Dirichlet boundary conditions
% Example 1 in Section 4 of 
% "Layer Potential Methods for Doubly-Periodic Harmonic Functions"
%
% This script uses the following classes/functions:
% - @segment, @quadr, @pointset, and selected utilities from @util in the 
%   mpspack library (by Alex Barnett and Timo Betcke, 2008–2012)
% - Custom methods written by Bohyun Kim (University of Utah, Math, 2025)

clear all; close all; clc;

verb = 0;
addpath('mps_helpers');  % Add helper functions to path

%% Define the geometry: one circular hole
M = 1;                          % Number of obstacles (here, one circle)
segs = cell(M, 1);             % Cell array for segment objects
cs = [0.5 + 0.5i];             % Center of the circle
rs = [0.1];                    % Radius
rphi = {@(t) -sin(t)};         % Dirichlet boundary data

% Define periodicity (tau = complex lattice vector)
tau = 1i;                     % Option: vertical periodicity
%tau = 0.5 + sqrt(3)/2 * 1i;    % Option: 60-degree lattice

%% Discretize the boundary
ss = segments.circles(cs, rs);          % Create circular segment
ss.requadrature([50]);                 % Coarse discretization
rhss = ss.copy;                         % Fine discretization for RHS eval
rhss.requadrature(250);

%% Build layer potential representation
ml = mylayerpot(ss, tau);                         % Operator using coarse grid
rhs = mylayerpot(rhss, tau, 'ms', rphi, 't');     % RHS using fine grid
rhseval = rhs.evalbd;                             % Evaluate RHS on boundary

% Map fine discretization to coarse for solving
conv_indxs = segments.coarse_to_fine_indxs(rhss, ss);

% Solve Dirichlet problem
pr = myproblems(ml, 'D');                  % 'D' = Dirichlet
pr.setuprhs(rhseval(conv_indxs));          % Set RHS
pr.solve;                                  % Solve linear system

%% Plot the solution in the fundamental domain
pt = Ptau(tau);  % Setup parallelogram domain

% Choose contour levels
levels = 10 * [-0.0035, -0.003, -0.0025, -0.002, -0.0015, -0.001, -0.0008, ...
               -0.0006, -0.0004, -0.0002, -0.0001, 0.0001, 0.0002, 0.0004, ...
                0.0006, 0.0008, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035];

% Plot solution with specified contour levels and black colormap
pr.plot(pt, levels, [0 0 0]);
