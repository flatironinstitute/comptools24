%WORKSHEET1_GEOMETRY_SOLVED
%
% NOTE: these are the solutions of the worksheet! It's recommended to try
% the worksheet without consulting these solutions first.
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - Representing curves using a panelization (charts defined on
%   intervals)
%   - The chunkIE approach to representing more complicated geometries as 
%   chunkgraph objects
%   - Defining integral kernels for use with chunkIE matrix builders and
%   post-processors
%   - Testing kernels with Gauss's identity and Green's identity
%

%% (Section 1) Legendre nodes 
%
% Legendre nodes are points in [-1,1] which are good for interpolating and
% integrating functions defined on [-1,1]
% 

%% 
% Legendre nodes and polynomials

close all; clearvars; clc;

% get some quantities for working with Legendre polynomials of order k

k = 16;
[x,w,u,v] = lege.exps(k);

% x -> Legendre nodes
% w -> Legendre weights
% u -> kxk matrix maps values at points x to Legendre series coefficients
%      of polynomial interpolant
% v -> kxk matrix evaluates Legendre series at points x. note v = inv(u)

figure(1)
clf
plot(x,zeros(size(x)),'rx')
hold on 

% can get Legendre polynomial values at various points from lege.pols

xviz = linspace(-1,1,200);
pviz = lege.pols(xviz(:),k);
plot(xviz,pviz(2,:))
plot(xviz,pviz(k/2,:))
plot(xviz,pviz(k+1,:))
legend('','P_1','P_{k/2}','P_{k}')

%%
% we can find the interpolant of f = sin(x) and look at the coefficients


% plot sin(n(x-0.1)^2) at 200 equispaced points in [-1,1] and plot values at
% Legendre nodes

n = 5;
ffun = @(t) sin(n*(t-0.1).^2);

fvals = ffun(x);
fviz = ffun(xviz);

figure(2)
clf
plot(xviz,fviz,'k-')
hold on
plot(x,fvals,'bo')

% compute values of interpolant at these points and compare to interpolant

fcoefs = u*fvals;
interpmat = pviz(1:k,:).';
finterp = interpmat*fcoefs;
plot(xviz,finterp,'b--')

% plot coefficients on a semilogy plot

figure(3)
clf
semilogy(abs(fcoefs),'b-o')

%%
% the integral obtained from sum(fvals.*w) is accurate for small n

fint_true = integral(ffun,-1,1,"RelTol",1e-15);
fint_approx = sum(fvals.*w);
abs(fint_true-fint_approx)

%% (Exercise 1)
%
% - use the Legendre nodes and weights to compute the integral of a 
% the function f = log(|t-1.1|) over [-1,1]
% - compare to an accurate integral obtained from MATLAB's integral
% function
% - see how the error depends on the order k 

ffun = @(t) log(abs(t-1.1));

fint_true = integral(ffun,-1,1,"RelTol",1e-15);

khigh = 30;
interrs = zeros(khigh,1);
for k = 1:khigh
    [x,w] = lege.exps(k);
    fvals = ffun(x);
    fint_approx = sum(fvals.*w);
    interrs(k) = abs(fint_approx-fint_true);
end

figure(4)
clf
semilogy(1:khigh,interrs,'b-o')

  

%% (Section 2) chunker geometries
%
% we represent a curve C by a collection of "chunks". Each chunk is
% represented as a Legendre series interpolant locally.

close all; clearvars; clc;

% we can build a chunker representation of a geometry given a function
% handle

fcurve = @(t) [cos(t(:).'); sin(t(:).')]; % should be 2 x numel(t) array

chnkr = chunkerfunc(fcurve);

figure(1)
clf
plot(chnkr)
hold on
plot(chnkr,'bo')
quiver(chnkr,'r')

%%
% curve positions for chunk j are stored in chnkr.r(:,:,j)
% derivatives are stored in chnkr.d(:,:,j) and chnkr.d2(:,:,j)
% normal vector is stored in chnkr.n(:,:,j)

% can look at a single panel
figure(2)
clf
plot(chnkr.r(1,:,1),chnkr.r(2,:,1),'rx')
hold on
quiver(chnkr.r(1,:,1),chnkr.r(2,:,1),chnkr.n(1,:,1),chnkr.n(2,:,1))

%%
% by default, uses k=16 nodes per panel, can change it up

pref = []; pref.k = 4;
chnkr2 = chunkerfunc(fcurve,[],pref);
chnkr2.nch
cparams = []; cparams.eps = 1e-2;
chnkr3 = chunkerfunc(fcurve,cparams,pref);
chnkr3.nch

%%
% can manipulate domains 

chnkr2 = chunkerfunc(@(t) starfish(t));

% move can do translations, rotations, etc

r0 = [0;0]; % shifts by -r0 before rotation and scaling
trotat = pi/4; % angle to rotate
scale = 1; % scales by scale
r1 = [3;4]; % shifts by +r1 after rotation and scaling
chnkr3 = reverse(chnkr2.move(r0,r1,trotat,scale));

figure(3)
clf
quiver(chnkr2,'r')
hold on
quiver(chnkr3,'b')

%%
% other functions/builders 

% chnk.curve.bymode position by fourier series for radius

modes = 0.2*randn(1,9);
modes(1) = 1.1*sum(abs(modes(2:end)));

chnkr = chunkerfunc(@(t) chnk.curves.bymode(t,modes));

figure(4)
clf
plot(chnkr,'b-o')

%%
% chunkerpoly - a rounded polygon with given vertices 

verts = chnk.demo.barbell();
chnkr2 = chunkerpoly(verts);

hold on
plot(chnkr2,'r-x')

%%
% chunkerfunc and chunkerpoly are adaptive by nature.
% chunkerfuncuni allows you to pick a number of chunks (evenly spaced in
% parameter space). convenient for convergence tests

nch = 5;
pref = []; pref.k = 5;
chnkr = chunkerfuncuni(@(t) starfish(t),nch,[],pref);

figure(5)
clf
plot(chnkr,'g-d')
hold on

%%
% another useful function is refine, which can oversample the grid (split
% each chunk in half)

refopts = []; refopts.nover = 1;
chnkr = chnkr.refine(refopts);
plot(chnkr,'b-o')
chnkr = chnkr.refine(refopts);
plot(chnkr,'r-x')

%%

% you can merge chunkers for multiply connected domains.
% it's a good idea to reverse orientation so that the normals
% are consistent

chnkr1 = chunkerfunc(@(t) chnk.curves.bymode(t,2));
chnkr2 = reverse(chunkerfunc(@(t) chnk.curves.bymode(t,1,[0.2,-0.3])));

chnkr = merge([chnkr1,chnkr2]);

figure(5)
plot(chnkr)
hold on
quiver(chnkr)

%%
% area will give you the area bounded by chunker if orientations are 
% right 

area_comp = area(chnkr);
area_true = 3*pi;
abs(area_comp-area_true)

%%
% chunkerinterior will tell you which points are inside/outside of a curve

x1 = linspace(-2,2,200); [xx,yy] = meshgrid(x1,x1);
pts = [xx(:).'; yy(:).'];

in = chunkerinterior(chnkr,pts);
in2 = chunkerinterior(chnkr,{x1,x1}); % faster option on meshgrids

figure(5)
clf
plot(chnkr)
hold on
scatter(pts(1,in2),pts(2,in2))

%% (Exercise 2) build a chunker 
% 
% come up with a curve parameterization for chunkerfunc or use chunkerpoly
% to define a domain, or both, and plot!
%
% advanced option: defining a multiply connected domain
%
% - combine the above routines to define a multiply connected domain
% consisting of a large, rounded rectangle with several smaller domains cut
% out of it
% - do a quiver plot of the domain to make sure you got the orientations
% right and that the boundaries don't intersect 
% - warning : you won't be able to merge if the chunkers have different 
% orders
% - set up a grid of points and find the points which are inside. make a
% plot of the inside points to check your work


W = 8; H = 8;
verts = [W/2 W/2 -W/2 -W/2; -H/2 H/2 H/2 -H/2];

chnkr0 = chunkerpoly(verts);

chnkr1 = reverse(chunkerfunc(@(t) starfish(t)));

scale = 0.3;

chnkrs = [];
for i = -3:3
    for j = -3:3
        % center on grid, perturb a little
        ctr = [i;j] + 0.1*scale*randn(2,1);
        % shift, scale, rotate by random amount
        chnkr2 = chnkr1.move([],ctr,2*pi*randn(),scale);
        chnkrs = [chnkrs,chnkr2];
    end
end

chnkrs = [chnkr0,chnkrs];
chnkr = merge(chnkrs);

figure(4)
clf
plot(chnkr)
hold on
quiver(chnkr)

x1 = linspace(-5,5,200); [xx,yy] = meshgrid(x1,x1);
pts = [xx(:).'; yy(:).'];

in2 = chunkerinterior(chnkr,{x1,x1});

figure(5)
clf
plot(chnkr)
hold on
scatter(pts(1,in2),pts(2,in2))


	 
%% (Exercise 3) if extra time: chunking ellipses
%
% - define the parameterization of an ellipse with semi-major axis a
% and semi-minor axis b. 
% - try chunking the ellipse with moderate order, say 8.
% - use chunkerfuncuni to produce chunkers with a range of nch values 
% - plot the error in the area of the chunker obtained and compare to 
% the expected error (depending on your implementation, you are likely to
% lose an order of accuracy because the chunker normal depends on
% derivatives of position)
% - loglog is a good plotting choice to see convergence order
%
% note: the area of an ellipse is pi*a*b if a and b are the semi-axis
% lengths
% 


a = 4;
b = 2;
true_area = pi*a*b;

ellipsefun = @(t) [a*cos(t(:).'); b*sin(t(:).')];

% a more advanced method to gain order of accuracy is to 
% also tell chunkerfunc about derivative info
ellipsefun2 = @(t) deal([a*cos(t(:).'); b*sin(t(:).')],...
    [-a*sin(t(:).'); b*cos(t(:).')]);

pref = []; pref.k = 8;

nchtotal = 10;
area_comps = zeros(nchtotal,1);
area_comps2 = zeros(nchtotal,1);
for nch = 1:nchtotal
    chnkr = chunkerfuncuni(ellipsefun,nch,[],pref);
    area_comps(nch) = area(chnkr);
    chnkr = chunkerfuncuni(ellipsefun2,nch,[],pref);
    area_comps2(nch) = area(chnkr); % these areas are surprisingly good!
end

nchs = 1:nchtotal;

figure(1)
clf
loglog(nchs,abs(area_comps-true_area),'b-x')
hold on
loglog(nchs,abs(area_comps2-true_area),'g-x')
loglog(nchs,nchs.^(-pref.k+1),'r--')    

  
%% (Section 3) more complicated geometries with chunkgraphs
%
% we describe problems with corners and multiple junctions using a
% "chunkgraph" structure. The chunkgraph consists of:
% - a collection of vertices, i.e. points where smooth curves meet. 
% Often, these are geometric singularities, like a corner or multiple 
% junction but can also be points on a smooth curve where mixed bounday
% conditions meet.
% - a collection of edges given by smooth curves. all vertices must be an 
% end point of an edge
%

close all; clearvars; clc;

% verts -> 2 x nverts array of vertex locations 
% edgesendverts => 2 x nedges array of connectivity information, each
% column indicates the start and end vertex for each edge 

% without further info, the constructor builds a polygonal domain with the
% given vertices and edges

verts = [1 0 -1 0; 0 1 0 -1]; edgesendverts = [1:3, 3, 4; 2:3, 1, 4, 1];
cg = chunkgraph(verts,edgesendverts);

% plot_regions will show the region ID numbers and edge numbers 

plot_regions(cg)
hold on
quiver(cg)

%% 
% you can specify curves for the edges, in place of straight lines 

nedge = 6;
thetas = linspace(0,2*pi,nedge+1); thetas = thetas(1:end-1);
verts = [cos(thetas(:).'); sin(thetas(:).')];
edgesendverts = [1:nedge; 2:nedge, 1];

fchnks = cell(nedge,1); 
cparams = []; cparams.ta = 0; cparams.tb = 2*pi;
for j = 1:nedge
    if mod(j,2)
        fchnks{j} = @(t) [t(:).'; sin(t(:).')];
    else
        fchnks{j} = [];
    end
end

% chunkgraph will automatically rotate, translate, scale to snap to
% vertices
cg = chunkgraph(verts,edgesendverts,fchnks,cparams);

figure(2)
clf
plot_regions(cg)
hold on
quiver(cg)

%% (Exercise 4) build a chunkgraph
%
% define your own chunkgraph with some straight and some curved sides 
% plot it to make sure you're getting what you imagined. 
%
