% Dirichlet Problem: Multiple Holes Example
% Solves for a basis spanning the nullspace of K*-I/2
% "Layer Potential Methods for Doubly-Periodic Harmonic Functions"
%
% This script uses the following classes/functions:
% - @segment, @quadr, @pointset, and selected utilities from @util in the 
%   mpspack library (by Alex Barnett and Timo Betcke, 2008–2012)
% - Custom methods written by Bohyun Kim (University of Utah, Math, 2025)

clear all; close all; clc; verb = 0; 
addpath("mps_helpers");

M = 2; segs = cell(M,1);
cs = [0.7+0.5*1i;0.3+ 0.3*1i]; rs = [0.1;0.15];
tau = 1i; %square torus example
Ns = [150; 150];

ss = segments.circles(cs,rs);
ss.requadrature(Ns);
ml = mylayerpot(ss,tau,'s');
ml.Dsoperator;
Atil = ml.Ds-0.5*eye(ml.ss.tN); % build matrix for K*-1/2
[V,D] = eig(Atil);
[eigenvals,eindx] = sort(diag(D)); 
pt = Ptau(tau);  % Setup parallelogram domain
xcount = 100; ycount = 100;
pt.add_grid(xcount,ycount);

% extracting eigenvector corresponding to zero eigenvalue
eigk = 1; 
eigkindx = eindx(eigk);
eigenval = eigenvals(eigk); 
psi= V(:,eigkindx);% get eigenvalue for zero.
ml.update_density(psi);
zz = pt.solgrid;
[~,seval] = ml.slayer(zz);

%plotting the colormap
h = pcolor(real(zz),imag(zz),10*seval); hold on;
set(h, 'EdgeColor', 'none');
colorbar;
hold on;
ss.addPtau(pt);
ss.plot2(1);
axis equal;
axis off;

%% plotting psi
figure;
for i = 1:M
    plot(2*pi*ss.ts(ss.indxs{i}),psi(ss.indxs{i}),'LineWidth',2); hold on;
end
t2 = xlabel("$\arg(z-a_i)$");
t1 = ylabel("$\psi_1$");
xlim([0,2*pi]);
set(t1,'Interpreter','latex');
set(t2,'Interpreter','latex');
ticksval = [0,pi/2,pi,3*pi/2,2*pi];
xticks(ticksval);
t3 = {"0","\pi/2","\pi","3\pi/2","2\pi"};
xticklabels(t3)

lg = legend('$D_1$','$D_2$');
set(lg,'Interpreter','latex');
ax1 = gca;
set(gca,'fontsize',17);
ylim([-0.1, 0.15])
