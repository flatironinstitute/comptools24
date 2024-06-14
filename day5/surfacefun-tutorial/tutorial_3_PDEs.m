%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partial differential equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Surfacefun package is capable of solving general variable-coefficient,
% second-order, linear, elliptic partial differential equations on open and
% closed surfaces. Consider such a PDE on a given surface \Gamma
%
%    L u(x) = f(x)
%
% where L is a an elliptic partial differential operator of the form
%
%  L u(x) = sum_{i=1}^3 sum_{j=1}^3 a_{ij}(x) ∂_i ∂_j u(x)
%         + sum_{i=1}^3 b_i(x) ∂_i u(x)
%         + c(x) u(x)
%
% with smooth coefficients a_{ij}, b_i, and c. Here, we identify ∂_1, ∂_2, ∂_3
% with the Cartesian x-, y-, and z-components of the  surface gradient,
% respectively. If Γ is not a closed surface, the PDE may also be subject to
% boundary conditions, e.g., u(x) = g(x) for x ∈ ∂Γ and some function g.

% Partial differential operators in Surfacefun are specified as MATLAB structs
% with properties defining the coefficients on each term appearing in L. The
% coefficients on the second-order terms are specified by setting:
% 
%   - pdo.dxx, pdo.dyy, pdo.dzz
%   - pdo.dxy, pdo.dyx
%   - pdo.dyz, pdo.dzy
%   - pdo.dxz, pdo.dzx
%
% Coefficients on the first-order terms can be set via:
%
%   - pdo.dx
%   - pdo.dy
%   - pdo.dz
% 
% The zeroth-order coefficient is specified via:
% 
%   - pdo.c
%
% Each coefficient may be a constant, a function handle, or a surfacefun. For
% instance, the Laplace-Beltrami operator can be specified via:

pdo = [];
pdo.dxx = 1;
pdo.dyy = 1;
pdo.dzz = 1;

% or more simply as:

pdo = [];
pdo.lap = 1;

% This sets the coefficients on the terms , , and to one and the rest to zero.

%% Let’s solve a simple Laplace–Beltrami problem on the surface of the sphere.
%  Since the spherical harmonics are eigenfunctions of the Laplace–Beltrami
%  operator on the sphere, we can construct a test problem analytically:

% Make a sphere mesh
p = 16;
nref = 2;
dom = surfacemesh.sphere(p+1, nref);

% Construct the true solution and corresponding righthand side
l = 3; m = 2;
sol = spherefun.sphharm(l, m);
sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
f = -l*(l+1)*sol;

% Specify the partial differential operator to be Laplace-Beltrami
pdo = [];
pdo.lap = 1;

%% To solve the PDE, we create a surfaceop. The surfaceop object encapsulates
%  the factorization and solution routines corresponding to the fast direct
%  solver.

L = surfaceop(dom, pdo, f)

%% It turns out that the Laplace–Beltrami problem on a closed surface is rank
%  deficient (by one). However, it is uniquely solvable under the mean-zero
%  conditions
%
%     integral(u) = integral(f) = 0
%
%  We can impose this condition by setting:

L.rankdef = true;

%% Now we can solve the PDE:

u = L.solve();
clf, plot(u), shg

%% Let's check the error:

norm(u - sol)

%% Surface PDEs on surfaces of arbitrary genus may be solved using surfaceop.
%  For example, here is the solution to a variable-coefficient surface Helmholtz
%  equation on a genus-1 stellarator geometry:

% Construct the stellarator geometry
p = 16; nu = 8; nv = 24;
dom = surfacemesh.stellarator(p+1, nu, nv);

% Variable-coefficient surface Helmholtz equation
pdo = [];
pdo.lap = 1;
pdo.c = @(x,y,z) 300*(1-z);

f = -1;
L = surfaceop(dom, pdo, f);
u = L.solve();

clf
plot(u), colorbar
shg

%% Now let's solve a problem on an open surface. We'll create an open surface by
%  extracting a subset of the patches from a closed surface:

rng(0)
p = 16;
nref = 2;
dom = surfacemesh.blob(p+1, nref);
dom = surfacemesh(dom.x(1:16), dom.y(1:16), dom.z(1:16));
clf
plot(dom), view(-110, 30), camlight
shg

%% We construct a surfaceop on an open surface in the same way as on a closed
%  surface, except now the L.solve() method requires Dirichlet boundary data to
%  be passed as an argument:

% Specify the righthand side and Dirichlet boundary data
f = -1;
bc = 0;

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
u = L.solve(bc);

clf
plot(u), view(-110, 30), colorbar
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modifying an existing surfaceop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A surfaceop is a direct solver for the specified surface PDE. This means that
% once a factorization of the problem is constructed, the factorization may be
% reused to solve the same PDE with different righthand sides and boundary data
% in a manner that is more efficient than creating a new surfaceop again and
% again.

%% To this end, a surfaceop may be factorized before it is given any data. This
%  is performed implicitly when L.solve() is called, but may be performed
%  explicitly by calling L.build():

L = surfaceop(dom, pdo);
L.build();

%% The surfaceop can now be given any righthand side and will efficiently update
%  its factorization accordingly:

L.rhs = @(x,y,z) sin(x.*y);
u = L.solve(bc);
clf
plot(u), view(-110, 30), colorbar
shg
