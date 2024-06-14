%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: Time-dependent reaction-diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Consider a nonlinear time-dependent PDE of the form
%
%    du/dt = Lu + N(u)
%
% with L a linear surface differential operator and N a nonlinear,
% non-differential operator. Such PDEs appear as models for reaction–diffusion
% systems, where L contains the diffusion terms and N the reaction terms. As the
% timescales of diffusion and reaction are often orders of magnitude different,
% explicit timestepping schemes can suffer from a severe time step restriction.
% Implicit or semi-implicit timestepping schemes can alleviate stability issues,
% e.g., by treating the diffusion term L implicitly.
%
% Fix a time step dt > 0. Discretizing in time with a simple backward Euler
% scheme results in a steady-state PDE at each time step of the form
%
%    (I - dt*L) u^{k+1} = u^k + dt*N(u^k)
%
% If the operator L is time-independent, then a solver for this PDE may be
% constructed once and reused at each time step with a different righthand side.
% We'll construct a surfaceop which does just that.

%% Complex Ginzburg-Landau equation

% The complex Ginzburg–Landau equation on a surface Γ is given by
%
%    du/dt = delta*(1+alpha*i)*lap(u) + u - (1+beta*i)*u*|u|^2
%
% where alpha, beta, and delta are parameters.

% Set parameters
dt = 0.03;
tend = 60;
alpha = 0;
beta = 1.5;
delta = 5e-4;

% Define nonlinear operator
N = @(u) u - (1+beta*1i)*u.*(abs(u).^2);

% Import the cow mesh
root = fileparts(fileparts(which('surfacefun')));
file = fullfile(root, 'models', 'cow.csv');
dom = surfacemesh.import(file, 'rhino');

% Construct linear solver
pdo = struct('lap', -dt*delta*(1+alpha*1i), 'c', 1);
L = surfaceop(dom, pdo);
L.build()

% Initial conditions
rng(1)
f = randnfun3(0.2, boundingbox(dom));
u = surfacefun(@(x,y,z) f(x,y,z), dom);

doplot = @(u) chain(@() clf, ...
                    @() plot(real(u)), ...
                    @() axis('off'), ...
                    @() material('dull'), ...
                    @() lighting('gouraud'), ...
                    @() colormap('turbo'), ...
                    @() colorbar, ...
                    @() clim('auto'), ...
                    @() view(180, -85), ...
                    @() camorbit(20, 0, 'data', 'y'), ...
                    @() camorbit(10, 0, 'data', 'x'), ...
                    @() camorbit(-3, 0, 'data', [0 0 1]), ...
                    @() camlight('headlight'));

doplot(u)
shg

% Simulate
t = 0;
kend = floor(tend/dt);
for k = 1:kend
    tic
    L.rhs = u + dt*N(u);
    u = L.solve();
    toc
    t = t + dt;
    if ( mod(k, 10) == 0 )
        fprintf('k = %d\n', k);
        doplot(real(u)), drawnow
    end
end

%% Turing system

% Consider the two-species reaction–diffusion system on a surface Γ given by,
%
%    du/dt = delta_u*lap(u) + alpha*u*(1-tau_1*v^2) + v*(1-tau_2*u),
%    dv/dt = delta_v*lap(v) + beta*v*(1+alpha*tau_1/beta*u*v) + u*(gamma+tau_2*v).
%
% Solutions u and v to this system can exhibit Turing patterns---namely, spots
% and stripes---depending on the choice of parameters delta_u, delta_v, alpha,
% beta, gamma, tau_1, and tau_2.

% Set parameters
dt = 0.1;
tend = 200;
delta_v = 1e-3;
delta_u = 0.516*delta_v;
alpha = 0.899;
beta = -0.91;
gamma = -0.899;
tau1 = 0.02;
tau2 = 0.2;

% Define nonlinear operators
Nu = @(u,v) alpha*u.*(1-tau1*v.^2) + v.*(1-tau2*u);
Nv = @(u,v) beta*v.*(1+alpha/beta*tau1*u.*v) + u.*(gamma+tau2*v);

% Import the cow mesh
root = fileparts(fileparts(which('surfacefun')));
file = fullfile(root, 'models', 'cow.csv');
dom = surfacemesh.import(file, 'rhino');

% Construct linear solvers
pdo = struct('lap', -dt*delta_u, 'b', 1);
Lu = surfaceop(dom, pdo);
Lu.build();
pdo = struct('lap', -dt*delta_v, 'b', 1);
Lv = surfaceop(dom, pdo);
Lv.build();

% Initial conditions
rng(1)
bb = boundingbox(dom);
f = randnfun3(0.2, bb);
u = surfacefun(@(x,y,z) f(x,y,z), dom);
f = randnfun3(0.2, bb);
v = surfacefun(@(x,y,z) f(x,y,z), dom);

doplot(u)
shg

% Simulate
t = 0;
kend = floor(tend/dt);
for k = 1:kend
    tic
    Lu.rhs = u + dt*Nu(u,v);
    Lv.rhs = v + dt*Nv(u,v);
    u = Lu.solve();
    v = Lv.solve();
    toc
    t = t + dt;
    if ( mod(k, 100) == 0 )
        fprintf('k = %d\n', k);
        doplot(real(u)), drawnow
    end
end
