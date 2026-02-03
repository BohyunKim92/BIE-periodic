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

ss = segments.neumann_regions;
Ns = 50*ones(ss.M,1);
ss.requadrature(Ns);
as = [1,2,3,4,-5,-2,-3,0];
tau = 0.5+(sqrt(3)/2)*1i;
ml = mylayerpot(ss,tau,'d');
ml2= mylayerpot(ss,tau,'ms');
pr = myproblems(ml,'N');
rhs = rhs1(ss,as,tau);
pr.setuprhs(rhs);
pr.solve;

pt = Ptau(tau); % set up parallelogram domain
xcount = 100; ycount = 100;
pt.add_grid(xcount,ycount)
levels = [linspace(-5,4,45),1e-6];
pr.plot(pt,levels)

function rhs = rhs1(ss,as,tau)
% test neumann bvp data
zshifts = repmat(ss.zs,[1,ss.M])- repmat(ss.cs.',[ss.tN, 1]);
[Gx,Gy] = myutils.gradG(zshifts,tau);
g = real(conj(repmat(ss.nus,[1,ss.M])).*(Gx+1i*Gy)); %dG/dnu(z)
rhs = g*as.';
end