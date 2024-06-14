%#ok<*NOPTS>
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before we can compute with functions and solve differential equations on
% surfaces, we have to define the surface itself. The surfacemesh class
% encapsulates the definition of the geometry.

% A surfacemesh consists of a collection of high-order quadrilateral patches
% whose union define the geometry. We choose to represent all patches of the
% surface with respect to their embedding space—the three-dimensional Cartesian
% coordinate system (x,y,z):

% Derived quantities such as metric tensor information and Jacobian factors are
% automatically computed and abstracted away from the user.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Surface representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each patch of a surfacemesh is represented using a set of high-order nodes—--
% specifically tensor-product Chebyshev nodes.

% Let's construct the set of tensor-product Chebyshev nodes in [-1,1]^2 and plot
% them:

n = 16;
[u, v] = chebpts2(n);
plot(u, v, 'ko', markerfacecolor='k')
axis equal off
shg

%% Now let’s create a surfacemesh with a single patch using these nodes. The
%  constructor takes as input the coordinates of the nodes on each patch as a
%  MATLAB cell array whose length is equal to the number of patches in the mesh:

x = u;
y = v;
z = cos(u).*cos(v);
dom = surfacemesh({x}, {y}, {z});
plot(dom), camlight
shg

%% We can verify that the mesh has a single element by asking for the length of
%  the mesh:

length(dom)

%% Now let's make a surface consisting of a few elements arranged in a 4 x 4
%  grid and map them to the graph of a given function.

mx = 4;
my = 4;

x = cell(mx*my, 1);
y = cell(mx*my, 1);
z = cell(mx*my, 1);

% Create a random function to define the surface
rng(0)
f = 0.2*randnfun2;

k = 1;
for i = 1:mx
    for j = 1:my
        % Map the Chebyshev nodes to smaller squares inside [-1,1]^2
        % and evaluate the given function
        uk = (u+1)/mx + 2/mx*(i-1) - 1;
        vk = (v+1)/my + 2/my*(j-1) - 1;
        x{k} = uk;
        y{k} = f(uk,vk);
        z{k} = vk;
        k = k+1;
    end
end

% Construct the surfacemesh
dom = surfacemesh(x, y, z);

plot(dom), camlight
shg

%% The mesh now has multiple patches:

length(dom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Built-in surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The surfacemesh class provides a number of built-in surface generation tools
% for creating simple surfaces. We'll start by creating a "cubed sphere" mesh,
% which consists of a cube mesh that has been inflated to live on the
% sphere.

%%  Create a sphere mesh of order p with two levels of refinement:

p = 16;
nref = 2;
dom = surfacemesh.sphere(p+1, nref);
plot(dom), camlight
shg

%% We can deform this mesh in various ways. For instance, we can create a
%  blob-like surface mesh by radially perturbed the nodes of the spherical mesh
%  according to a random smooth function. This is encapsulated in the
%  surfacemesh.blob routine.

% Create a blob mesh of order p with two levels of refinement:
rng(0)
dom = surfacemesh.blob(p+1, nref);
plot(dom), camlight
shg

%% If we keep calling this function, we'll generate new surfaces due to the
%  randomness of the algorithm:

rng(0)
for k = 1:3
    dom = surfacemesh.blob(p+1, nref);
    figure(k), clf, plot(dom), camlight
end
alignfigs

%% Surface meshes of any genus are supported. Here is a smoothly deformed torus:

p = 16; nu = 8; nv = 24;
dom = surfacemesh.torus(p+1, nu, nv);
plot(dom), camlight
shg

%% The mesh does not even have to be smooth between patches...

p = 16; nu = 4; nv = 32;
dom = surfacemesh.twisted_torus(p+1, nu, nv);
plot(dom), camlight
shg

%% ...or even orientable!

p = 16; nu = 30; nv = 7;
dom = surfacemesh.mobius(p+1, nu, nv);
plot(dom), camlight
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Importing an existing mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CAD packages such as Gmsh or Rhino can be used to create a high-order smooth
% quadrilateral mesh. Surfacefun supports importing meshes from these packages:

root = fileparts(fileparts(which('surfacefun')));
file = fullfile(root, 'models', 'cow.csv');
dom = surfacemesh.import(file, 'rhino');
plot(dom)
view(180, -85)
camorbit(20, 0, 'data', 'y')
camorbit(10, 0, 'data', 'x')
camorbit(-3, 0, 'data', 'z')
camlight
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convenient surfacemesh routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 16;
nref = 2;
dom = surfacemesh.sphere(p+1, nref);

%% Only plot the wireframe of the surfacemesh, not the underlying surface:

plot(dom, surface='off')
shg

%% Make a mesh plot of the high-order nodes on each patch:

mesh(dom)
shg

%% The number of patches in the surfacemesh:

length(dom)

%% The polynomial order of each patch:

order(dom)

%% The total number of degrees of freedom in the surfacemesh:

numel(dom)

%% The volume enclosed by a closed surfacemesh:

volume(dom) - 4/3*pi

%% The surface area of the surfacemesh:

surfacearea(dom) - 4*pi

%% The smallest box that contains the surfacemesh:

boundingbox(dom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modifying a mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A mesh may modified by changing the polynomial degree used to represent each
% patch ("p-refinement") or by the changing the number of patches
% ("h-refinement").

%% Change the order of mesh by resampling to a new set of nodes. Let's resample
%  the mesh to use 3 nodes per patch:

dom2 = resample(dom, 3);
plot(dom2), camlight
shg

%% Now refine the mesh by dividing each patch into four:

dom3 = refine(dom2);
plot(dom3), camlight
shg
