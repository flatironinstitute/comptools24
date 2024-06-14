%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions and vector fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Once a surfacemesh has been constructed, we may define scalar functions and
% vector fields on the surface.

rng(0)
p = 16; nref = 2;
dom = surfacemesh.blob(p+1, nref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scalar functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The fundamental object which represents a scalar function on a surface is a
%  surfacefun. A surfacefun may be constructed from a function handle
%  representing a given function in Cartesian (x,y,z) coordinates on a given
%  surfacemesh:

f = surfacefun(@(x,y,z) cos(6*x).*y + exp(z), dom)

%% Let's plot the function:

clf
plot(f), hold on, plot(dom), colorbar
shg

%% Many standard MATLAB arithmetic functions have been overloaded:

x = surfacefun(@(x,y,z) x, dom);
g = abs(f + 2*x);
clf
plot(g), colorbar
shg

%% We can also visualize a surfacefun using a contour plot:

clf
contour(f, linewidth=2)
axis off
shg

%% We may numerically differentiate a function using the built-in diff or grad
%  routines, which automatically take into account the on-surface metric. For
%  example:

[fx, fy, fz] = grad(f);
plot([fx fy fz])

%% Higher-order derivatives may be constructed by composing these operations.
%  For example, here is the surface Laplacian---or the Laplace-Beltrami operator
%  ---applied to our function:

clf
plot(lap(f)), colorbar
shg

%% The definite integral of a function over the surface is given by:

integral(f)

%% Similarly, the mean of the function is the integral of the function divided
%  by the surface area:

mean(f)

%% We can estimate the minimum and maximum values of a function:

[minEst(f) maxEst(f)]

%% We can take the inner product of two surfacefuns as if they were column
%  vectors:

f' * g

%% Let's make some random functions. We will use Chebfun to construct random
%  smooth functions in 3D and then construct surfacefuns for the restriction of
%  these functions to the surface:

box = boundingbox(dom);
h = [];
for k = 1:6
    r = randnfun3(0.5, box);
    h = [h surfacefun(@(x,y,z) r(x,y,z), dom)];
end
plot(h)

%% They are not orthogonal:

h' * h(:,1)

%% ...but we can orthogonalize them:

H = orth(h);
H' * H
plot(H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Norms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The L^2 norm of a surfacefun may be computed via:

norm(f)

%% Other norms are implemented as well. The L^inf norm is computed via:

norm(f, inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vector fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The surfacefunv object represents a three-component vector field over a
% surfacemesh. Each component is itself represented as a scalar surfacefun.
% Let's make quiver plot of the normal vectors over our surface. We'll plot 6
% vectors per patch and scale their lengths by 0.2:

v = normal(dom);
clf
quiver(v, 0.2, 6)
shg

%% The surface gradient of a surfacefun is a surfacefunv:

grad(f)

%% The gradient is tangent to the surface, as we can see from a quiver plot:

clf
quiver(grad(f), 0.05, 6)
shg

%% The surface divergence of the surface gradient is equal to the surface
%  Laplacian:

norm(div(grad(f)) - lap(f))

%% The mean curvature of a surface can be related to its the normal vector field
%  via the surface divergence:

clf
kappa = div(v) / 2;
plot(kappa)
colorbar
shg

%% We can also take the surface curl of a surfacefunv:

v = surfacefunv(@(x,y,z) cos(2*x), ...
                @(x,y,z) sin(4*y), ...
                @(x,y,z) sin(3*z), dom);
clf
quiver(curl(v), 0.1, 6)
shg

%% Cross products are also supported:

clf
w = cross(curl(v), [0 1 0]);
plot(norm(w)), hold on
quiver(w, 0.1, 8, 'k'), hold off
colorbar
shg

%% ...as well as dot products:

clf
plot(dot(curl(v), [0 1 0]))
shg
