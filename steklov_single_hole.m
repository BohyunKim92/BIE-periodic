% Steklov Problem: Single Hole Example
% Solves Steklov Eigenvalue problem
% Example 4 in Section 4 of 
% "Layer Potential Methods for Doubly-Periodic Harmonic Functions"
%
% This script uses the following classes/functions:
% - @segment, @quadr, @pointset, and selected utilities from @util in the 
%   mpspack library (by Alex Barnett and Timo Betcke, 2008–2012)
% - Custom methods written by Bohyun Kim (University of Utah, Math, 2025)

clear all; close all; clc; 
addpath('mps_helpers');
%% Define the geometry: one circular hole
M = 1;                          % Number of obstacles (here, one circle)
cs = [0.5 + 0.5i];             % Center of the circle
rs = [0.2];                    % Radius

% Define periodicity (tau = complex lattice vector)
%tau = 1i;                     % Option: square torus
tau = 0.5 + sqrt(3)/2 * 1i;    % Option: equilateral torus

%% Discretize the boundary
ss = segments.circles(cs, rs);          % Create circular segment
ss.requadrature([50]);                 % Coarse discretization

ml = mylayerpot(ss,tau,'s');
pr = myproblems(ml, 'S');                  % 'S' = steklovEVP
pr.solve;


% calculate errors
[eigenvals,eindx] = sort(pr.evals); 
real_eig = sort(real(eigenvals));
half_real_eig = 0.5*real_eig; % to consider scaling
if tau == 1i
    exact = [0;3.21737540790552735473880286001400036767774798208487;3.21737540790552735473880286001400036767774798208487;4.85099530552467697892257589130439715581461931719259;5.15358084940676223549771471754234765157435969419525;7.50305008416767542642635086056165243882709526430554;7.50305008416767542642635086056165243882709526430554];
elseif tau ==0.5 + sqrt(3)/2 * 1i
    exact = [0;3.34865594380260534169550288243470971962587318064277;3.34865594380260534169550288243470971962587318064277; 4.99978881548382813234141616969113198885117552416465; 4.99978881548382813234141616969113198885117552416465; 7.44392530690947308002824485738760008901145380307620; 7.55649710043624518482844840631875099119732734059433];
else
    error("eigenvalues comparable only for tau = 1i or tau = 0.5 + sqrt(3)/2 * 1i case.");
end

err = abs(half_real_eig(1:7)-exact);
imag_eig = imag(eigenvals);
if norm(imag_eig,"inf") < 10^(-11)
    disp("imaginary eigenvalue sufficiently small");
end
