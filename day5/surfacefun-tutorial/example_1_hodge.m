%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: Hodge decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Laplace–Beltrami problem arises when computing the Hodge decomposition of
% tangential vector fields. For a vector field f tangent to the surface Γ, the
% Hodge decomposition writes f as the sum of curl-free, divergence-free, and
% harmonic components,
%
%    f = grad(u) + n x grad(v) + w
%
% where u and v are scalar functions on Γ and w is a harmonic vector field,
% i.e.,
%
%    div(w) = 0,   div(n x w) = 0.
%
% Such vector fields play an important role in integral-equation-based methods
% for computational electromagnetics. To numerically compute such a
% decomposition, one may solve two Laplace–Beltrami problems for u and v,
%
%    lap(u) = div(f),   lap(v) = -div(n x f),
%
% and then set w = f - grad(u) - n x grad(v).


% Construct a toroidal mesh
p = 16; nu = 16; nv = 48;
dom = surfacemesh.torus(p+1, nu, nv);

% Make a random smooth tangential vector field
rng(0)
gx = randnfun3(10, boundingbox(dom));
gy = randnfun3(10, boundingbox(dom));
gz = randnfun3(10, boundingbox(dom));
g = cross([0 1 1], surfacefunv(@(x,y,z) gx(x,y,z), ...
                               @(x,y,z) gy(x,y,z), ...
                               @(x,y,z) gz(x,y,z), dom));
vn = normal(dom);
f = -cross(vn, vn, g);

% Compute the Hodge decomposition
tic
[u, v, w] = hodge(f);
toc

%% Plot the resulting fields

doplot = @(f,k,t) chain(@() figure(k), ...
                        @() clf, ...
                        @() plot(norm(f)), ...
                        @() hold('on'), ...
                        @() quiver(normalize(f), 0.6, 3, 'k'), ...
                        @() axis('equal', 'off'), ...
                        @() colormap('turbo'), ...
                        @() colorbar, ...
                        @() title(t, 'FontSize', 16));

doplot(f, 1, 'Vector field'), hold on, plot(dom)
doplot(grad(u), 2, 'Curl-free component')
doplot(cross(vn, grad(v)), 3, 'Divergence-free component')
doplot(w, 4, 'Harmonic component')
alignfigs

%% Let's check how numerically harmonic the resulting w field is:

norm(div(w))
norm(div(cross(vn, w)))
