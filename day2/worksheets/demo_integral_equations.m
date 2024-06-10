% Worksheet 2: Layer potential evaluators and solvers
%
%  Learning Goals:
%    - Compute hard(er) integral using layer potential evaluators
%    - Solve integral equations
%    - Verify order of accuracy in analytic test/self-convergence tests

%% Section 1: Layer potential evaluator
clear all
clear classes

%% Run startup script
run ~/git/fmm3dbie/matlab/startup.m %Enter your path to fmm3dbie startup ehre

% Redo Sprime but with integral equation evaluator
S = geometries.sphere(1, 2);

targinfo = [];
targinfo.r = [0.995; 0.01; 0.03];

densities = ones(S.npts,1);
p = -eval_fields(S, 'l', 'd', densities, targinfo, 1e-6);

fprintf('Error in potential=%d\n', abs(p-1)/sqrt(4*pi));

%% Exercise 1: Repeat exercise on ellipsoid
% Find close points inside and outside the ellipsoid and compute the
% laplace double layer potential for those points with constant densities

%% Section 2: Solve exterior Dirichlet boundary value problem

S = geometries.sphere(1, 2);
% Get boundary data
xyz_in = [0.3;0.5;0.1];  % Interior point
xyz_out = [1.3;-5.2;0.1]; % exterior point


zk = 1.1; % wave number
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk, src_info, S, 's');

rep_pars = [-1j*zk, 1];
eps = 1e-7;
sig = solver(S, 'helm', 'dir', rhs, eps, zk, rep_pars);


%% Section 3: Test accuracy of solution

targ_info = [];
targ_info.r = xyz_out;
p = eval_fields(S, 'helm', 'dir', sig, targ_info, eps, zk, rep_pars);

pot_ex = helm3d.kern(zk, src_info, targ_info, 's');
fprintf('Error in iterative solver=%d\n',abs(p-pot_ex)/abs(pot_ex));

%% Section 4: Explore all available solver options
help solver
%%
help eval_fields
%% Exercise 2: Estimate the order of convergence
% Do a convergence study


%% Exercise 3: Merge + verify order of convergence
S = geometries.sphere(1, 2);
S2 = translate(S, [3;0;0]);
S = merge([S, S2]);

close('all')
plot(S)


%% Exercise 4: repeat with favorite boundary value problem, 
% can use surface smoother geometries as well, convergence study
% + post processing routines

