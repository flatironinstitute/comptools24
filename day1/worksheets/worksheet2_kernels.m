%WORKSHEET2_KERNELS
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - Defining an integral kernel to work with curve data
%   - Visualizing kernels
%   - Checking kernel identities
%   

%% (Section 1) Integral kernels and curve data 
%
% A chunker object stores position information (r), derivative information
% (d), second derivative information (d2), and normal vectors (n) for each
% point of the discretization. This same information is interpolated and
% passed to integral kernels to do things like compute integrals, etc. 
% The derivatives are with respect to the panel discretization, not arc
% length.
%
% A kernel evaluator in chunkIE can expect two inputs (s,t) which have the 
% fields 
% s.r, t.r -> positions (size 2 x ns and 2 x nt arrays, respectively)
% s.d, t.d -> first derivatives (if defined)
% s.d2, t.d2 -> second derivatives (if defined)
% s.n, t.n -> normal vectors (if defined)
%
% It should return the matrix of interaction, a nt x ns matrix 

close all; clearvars; clc;
set(0, 'DefaultLineLineWidth', 2);

% see below for definition of singlelayer_eval
s = []; s.r = randn(2,2); t = []; t.r = randn(2,3);
mat = singlelayer_eval(s,t)

%% 
% with this definition, we can also plug in a chunker object as source or
% target

s = []; s.r = [1.5;1.5]; 
chnkr = chunkerfunc(@(t) starfish(t));

mat = singlelayer_eval(s,chnkr);

figure(1)
clf
plot(s.r(1),s.r(2),'rx')
hold on
plot(chnkr)

arc = arclengthfun(chnkr);

figure(2)
plot(arc(:),mat(:))

%% 
% there are several built in kernels in chunkIE, available as 
% kernel class objects

kerns = kernel('laplace','s') % has evaluator, singularity info, fmm info

%%

% for large numbers of sources/targets can use fmm to evaluate the quantity 
% mat*val where val is a vector of source strengths 

ns = 1000; nt = 1000;
s = []; 
s.r = randn(2,ns);
t = [];
t.r = randn(2,nt);

strength = randn(ns,1);

tic; mat = kerns.eval(s,t); % forming matrix densely 
u1 = mat*strength; toc % product directly 

eps = 1e-12; % precision for fmm
tic; u2 = kerns.fmm(eps,s,t,strength); toc % this calling sequence is standard

norm(u1-u2)

%%
%

% we can visualize by placing targets on a grid 

x1 = linspace(-2,2,200); [xx,yy] = meshgrid(x1,x1);
t = [];
t.r = [xx(:).'; yy(:).'];

% some random "dipole" charges 

ns = 4;
s = []; s.r = randn(2,ns); s.n = randn(2,ns);
strength = randn(ns,1);

zk = 10+1i;
kerndhelm = kernel('h','d',zk);

eps = 1e-3;
u = kerndhelm.fmm(eps,s,t,strength);

figure(3)
clf
h = pcolor(reshape(real(u),size(xx))); set(h,'EdgeColor','none');
clim([-3,3])
colormap(redblue)

%%
%

% you can see some of the most common available kernels in the docs 

help kernel

%% (Exercise 1)
%
% pick a built-in kernel and visualize it on a grid (note: fmm is not 
% defined for all kernels). 
% Warning: Some kernels are vector/matrix (see opdims
% array, if [1 1], then scalar, otherwise the kernel is vector/matrix
% valued).
%
% Some scalar options: Helmholtz single and double layer, Laplace single
% and double.

%% (Section 2) Integrating a kernel over a curve
%

% the chunker object has a "native" quadrature rule based on the Legendre
% integration rule on each panel (and the arc length density ds/dt =
% sqrt(x_1'^2 + x_2'^2)) 

chnkr = chunkerfunc(@(t) [cos(t(:).'); sin(t(:).')]);
wts = chnkr.wts;

% the sum of the weights gives the length of the curve 

abs(sum(wts(:)) - 2*pi)

%%
% we can integrate a function over the curve by sampling the function at
% the nodes, multiplying by the weights, and summing the result

% set up exponential functions on the circle
fcurve = @(r,n) (r(1,:) + 1i*r(2,:)).^n;

ns = 0:5;
fcurve_vals = zeros(chnkr.npt,length(ns));
for j = 1:length(ns)
    n = ns(j);
    fcurve_vals(:,j) = fcurve(chnkr.r,n);
end

% can get the integrals by multiplying by the weights

fints = fcurve_vals.'*chnkr.wts(:)

%%
% one way of testing kernels is to make sure they satisfy 
% certain identities 
%
% A version of Green's identity holds for several kernels
%
% u = S[du_b/dn] - D[u_b] where u_b is the boundary data for 
% some potential u. 
%
% The value of S[mu] is an integral over the boundary 
%
% S[mu](x) = int_Gamma G(x,y) mu(y) ds(y)
%

chnkr = chunkerfunc(@(t) starfish(t));

% define u as the field induced by some random sources *outside* the domain
ns = 5;
s = []; 
rads = 2 + 0.5*rand(1,ns); thetas = 2*pi*rand(1,ns);
s.r = [cos(thetas); sin(thetas)].*rads;
strength = randn(ns,1);

figure(4)
clf
plot(chnkr)
hold on
plot(s.r(1,:),s.r(2,:),'rx')

% get boundary data associated to these sources 

zk = 5;
kerns = kernel('h','s',zk);
kernsprime = kernel('h','sp',zk);
kernd = kernel('h','d',zk);

ubdry = kerns.eval(s,chnkr)*strength;
unbdry = kernsprime.eval(s,chnkr)*strength;

% then evaluate the integral at a target inside the curve using the 
% native integral rule 

t = []; t.r = -0.25+0.5*rand(2,1);
plot(t.r(1),t.r(2),'gx')

ucomp = kerns.eval(chnkr,t)*(chnkr.wts(:).*unbdry) - ...
    kernd.eval(chnkr,t)*(chnkr.wts(:).*ubdry);

utrue = kerns.eval(s,t)*strength;

abs(ucomp-utrue)

%% (Exercise 2) visualize error in Green's identity at many targets
% 
% - set up a grid of targets for visualizing the error
% - make sure sources are outside (can use the code above)
% - compare true formula for the potential to the Green's identity 
% - set the error to nan for targets which are not inside (use
% chunkerinterior to find targets inside/outside) for plotting purposes
% - where is the error highest?


%% (Section 3) corrected quadrature
%

% corrected quadrature is available via chunkerkerneval
% if defined for the kernel, chunkerkerneval will use the fmm to accelerate

chnkr = chunkerfunc(@(t) starfish(t),cparams);

% here, we'll check Gauss' ID using quadarature
% D[1] = -1 inside 

x1 = linspace(-2,2,200); [xx,yy] = meshgrid(x1,x1);
t = [];
t.r = [xx(:).'; yy(:).'];

in = chunkerinterior(chnkr,{x1,x1});

kernd = kernel('l','d');

dens = ones(chnkr.npt,1);

ucomp = chunkerkerneval(chnkr,kernd,dens,t.r(:,in(:)));
ucompsmooth = chunkerkerneval(chnkr,kernd,dens,t.r(:,in(:)),struct("forcesmooth",true));

err_plot= nan(size(xx));
errsmooth_plot = nan(size(xx));

err_plot(in) = abs(ucomp+1);
errsmooth_plot(in) = abs(ucompsmooth+1);

figure(5)
clf
tiledlayout(1,2,'TileSpacing','compact')
nexttile
h = pcolor(reshape(log10(err_plot),size(xx))); set(h,'EdgeColor','none');
colorbar
nexttile
h = pcolor(reshape(log10(errsmooth_plot),size(xx))); set(h,'EdgeColor','none');
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  define functions below here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = singlelayer_eval(s,t)
%SINGLELAYER_EVAL get the interaction matrix for a single layer kernel
%

% find dimensions 
ns = size(s.r(:,:),2);
nt = size(t.r(:,:),2);

xs = s.r(1,:); ys = s.r(2,:);
xt = t.r(1,:).'; yt = t.r(2,:).';

rr = (xs-xt).^2 + (ys-yt).^2;

mat = -log(rr)/(4*pi);
end