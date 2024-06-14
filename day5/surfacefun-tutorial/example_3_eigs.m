%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: Eigenvalue problems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Surfaceop can be used to solve eigenvalue problems. In fact, the EIGS()
% routine is overloaded to act on a surfaceop. Internally, the fast direct
% solver for a surfaceop is computed once and efficiently reused.

%% Eigenvalues on the sphere

% Create a surfaceop for the surface Laplacian on a sphere:
n = 16; nref = 2;
dom = surfacemesh.sphere(n, nref);

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo);
L.rankdef = true;

% Compute the 6 smoothest eigenfunctions and their corresponding eigenvalues:
tic
[V, D] = eigs(L);
toc

diag(D)
plot(V)

%% Eigenvalues on the Mobius strip

% Create a surfaceop for the surface Laplacian on a Mobius strip:
n = 16; nu = 30; nv = 6;
dom = surfacemesh.mobius(n, nu, nv);
plot(dom)

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo);

% Compute the 6 smoothest eigenfunctions and their corresponding eigenvalues:
tic
[V, D] = eigs(L);
toc

diag(D)
plot(V)

%% Find the 10 eigenfunctions with eigenvalues closest to -1000:

tic
[V, D] = eigs(L, 10, -1000);
toc

diag(D)
plot(V)
