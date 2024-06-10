%WORKSHEET3_STOKES
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - This worksheet demonstrates solving the Stokes equation using the
%   chunkIE package 
%   - guides the user through creating geometries, selecting an appropriate
%   Stokes kernel, building the system matrix for the BIE, solves the
%   system matrix for the BIE, and then plots the result
%
%   - Defining Stokes integral kernel to work with curve data
%   - Checking Stokes kernel identities
%   - Demonstrating a convergence test 
%   - Also demonstrating a flow problem in a pipe with obstacles inside 


%% (Section 1) Stokes integral kernels and curve data 
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
% It should return the matrix of interaction, a (2xnt) x (2xns) velocity
% matrix, and a nt x (2xns) pressure matrix

close all; clearvars; clc;
set(0, 'DefaultLineLineWidth', 2);

% see below for definition of interleaving stokes singlelayer_eval
mu = 1;
s = []; s.r = randn(2,2); t = []; t.r = randn(2,3);
matsvel = stk_singlelayer_eval(s,t,mu)
matspres = stk_singlelayer_pres_eval(s,t,mu)

%% 
% with this definition, we can also plug in a chunker object as source or
% target

s = []; s.r = [1.5;1.5]; 
chnkr = chunkerfunc(@(t) starfish(t));

matsvel = stk_singlelayer_eval(s,chnkr,mu);
matspres = stk_singlelayer_pres_eval(s,chnkr,mu);
sigma = [1/2; 1/2];

figure(1)
clf
plot(s.r(1),s.r(2),'rx')
hold on
plot(chnkr)

arc = arclengthfun(chnkr);

figure(2)
clf
tiledlayout(1,3,'TileSpacing','compact')
nexttile
plot(arc(:),matsvel(1:2:end,:)*sigma)
nexttile
plot(arc(:),matsvel(2:2:end,:)*sigma)
nexttile
plot(arc(:),matspres*sigma)


%% 
% chunkIE has built-in Stokes kernels (fmm coming soon!)

kerns = kernel('stok','svel',mu); % has evaluator, singularity info
kernspres = kernel('stok','spres',mu);

%%
% can use mat*val where val is a vector of source strengths 

ns = 1000; nt = 1000;
s = []; 
s.r = randn(2,ns);
t = [];
t.r = randn(2,nt);

strength = randn(ns,2);

tic; matsvel = kerns.eval(s,t); % forming matrix densely 
u1 = matsvel*reshape(strength',[],1); toc % product directly 

matsvel2 = stk_singlelayer_eval(s,t,mu);
u2 = matsvel2*reshape(strength',[],1); toc

norm(u1-u2)

tic; matspres = kernspres.eval(s,t); % forming matrix densely 
u1pres = matspres*reshape(strength',[],1); toc % product directly 

matspres2 = stk_singlelayer_pres_eval(s,t,mu);
u2pres = matspres2*reshape(strength',[],1); toc

norm(u1pres-u2pres)

%%
%

% we can visualize by placing targets on a grid 

xtar = linspace(-2,2,200); [xx,yy] = meshgrid(xtar,xtar);
t = [];
t.r = [xx(:).'; yy(:).'];

% some random Stokes doublet

ns = 4;
s = []; s.r = randn(2,ns); s.n = randn(2,ns);
strength = randn(ns,2);

kernd = kernel('stok','d',mu);
kerndpres = kernel('stok','dpres',mu);

matdvel = kernd.eval(s,t);
u = matdvel*reshape(strength',[],1);
matdpres = kerndpres.eval(s,t);
upres = matdpres*reshape(strength',[],1);

figure(3)
clf
tiledlayout(2,2,'TileSpacing','compact')
nexttile
h = pcolor(xx,yy,reshape(u(1:2:end),size(xx))); set(h,'EdgeColor','none');
clim([-3,3])
colormap(redblue)
nexttile
h = pcolor(xx,yy,reshape(u(2:2:end),size(xx))); set(h,'EdgeColor','none');
clim([-3,3])
colormap(redblue)
nexttile
streamslice(xx,yy,reshape(u(1:2:end),size(xx)),reshape(u(2:2:end),size(xx)))
xlim([-2,2]),ylim([-2,2])
nexttile
h = pcolor(xx,yy,reshape(upres,size(xx))); set(h,'EdgeColor','none');
clim([-3,3])
colormap(redblue)

%% (Section 2) Visualize error in Green's identity
%
% one way of testing kernels is to make sure they satisfy 
% certain identities 
%
% A version of Green's identity holds for several kernels
%
% u = S[trac_b] - D[u_b] where u_b is the boundary velocity and trac_b is 
% the boundary traction for some velocity field u. 
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
strength = randn(ns,2);

figure(1)
clf
plot(chnkr)
hold on
plot(s.r(1,:),s.r(2,:),'rx')

% get boundary data associated to these sources 

kerns = kernel('stok','s',mu);
kernstrac = kernel('stok','strac',mu);
kernd = kernel('stok','d',mu);

ubdry = kerns.eval(s,chnkr)*reshape(strength',[],1);
tracbdry = kernstrac.eval(s,chnkr)*reshape(strength',[],1);

% then evaluate the integral at a target inside the curve using the 
% native integral rule 

xtar = linspace(-2,2,300); [xx,yy] = meshgrid(xtar,xtar);
t = [];
t.r = [xx(:).'; yy(:).'];

in = chunkerinterior(chnkr,t);

ucomp = kerns.eval(chnkr,t)*(reshape([chnkr.wts(:),chnkr.wts(:)]',[],1).*tracbdry) - ...
    kernd.eval(chnkr,t)*(reshape([chnkr.wts(:),chnkr.wts(:)]',[],1).*ubdry);

utrue = kerns.eval(s,t)*reshape(strength',[],1);

errux = abs(ucomp(1:2:end)-utrue(1:2:end));
erruy = abs(ucomp(2:2:end)-utrue(2:2:end));

errplotx = nan(size(xx));
errplotx(in) = errux(in);
errploty = nan(size(xx));
errploty(in) = erruy(in);

figure(2)
clf
tiledlayout(1,2,'TileSpacing','compact')
nexttile
h = pcolor(xx,yy,log10(errplotx)); set(h,'EdgeColor','none');
colorbar
hold on
plot(chnkr,'k-o')
nexttile
h = pcolor(xx,yy,log10(errploty)); set(h,'EdgeColor','none');
colorbar
hold on
plot(chnkr,'k-o')


%% (Exercise 2) corrected quadrature
% corrected quadrature is available via chunkerkerneval
% - use the above grid of targets for visualizing the error
% - make sure sources are outside (can use the code above)
% - compare true formula for the potential to the Green's identity 
% - set the error to nan for targets which are not inside (use
% chunkerinterior to find targets inside/outside) for plotting purposes

% here, we'll check Green's identity using quadarature
% 

ucomp = chunkerkerneval(chnkr,kerns,tracbdry,t.r(:,in(:))) - ...
      chunkerkerneval(chnkr,kernd,ubdry,t.r(:,in(:)));

utrue = kerns.eval(s,struct('r',t.r(:,in(:))))*reshape(strength',[],1);


errux = abs(ucomp(1:2:end)-utrue(1:2:end));
erruy = abs(ucomp(2:2:end)-utrue(2:2:end));

errplotx = nan(size(xx));
errplotx(in) = errux;
errploty = nan(size(xx));
errploty(in) = erruy;

figure(3)
clf
tiledlayout(1,2,'TileSpacing','compact')
nexttile
h = pcolor(reshape(log10(errplotx),size(xx))); set(h,'EdgeColor','none');
colorbar
nexttile
h = pcolor(reshape(log10(errploty),size(xx))); set(h,'EdgeColor','none');
colorbar

%% (Section 3) Interior Stokes Dirichlet problem
%
% with the representation u = c*S[sigma]+D[sigma], the boundary integral
% equation is -1/2 sigma + c*S[sigma] + D[sigma] + W[sigma] = f, where f 
% is the Dirichlet boundary data. The operator W removes the nullspace 
%
% It has the value 
%
% W[sigma](x) = n(x) \int_{\Gamma} n(y) \cdot sigma(t) ds(y)

k = 2;
nch = 5*k;
chnkr = chunkerfuncuni(@(t) starfish(t),nch);

% define the combined layer Stokes representation kernel
c = -1;
mu = 1;
kerndvel = kernel('stok','dvel',mu);
kernsvel = kernel('stok','svel',mu);

kerncvel = kerndvel + c*kernsvel;

% get a matrix discretization of the boundary integral equation 
cmat = chunkermat(chnkr,kerncvel);
sysmat = cmat - 1/2*eye(2*chnkr.npt) + normonesmat(chnkr);

% we can define u as the field induced by some random sources *outside* the
% domain and get boundary data associated to these sources
ns = 5;
s = []; 
rads = 2 + 0.5*rand(1,ns); thetas = 2*pi*rand(1,ns);
s.r = [cos(thetas); sin(thetas)].*rads;
strength = randn(ns,2);
rhs = kerns.eval(s,chnkr)*reshape(strength',[],1);

% and solve
dens = sysmat\rhs;

%% Now plot the solution

% grid for plotting purposes
xts = linspace(-1.5,1.5);
yts = xts;
[xtar,ytar] = meshgrid(xts,yts);
sz = size(xtar);

int_inds = chunkerinterior(chnkr,{xts,yts});

targs = [];
targs.r = [xtar(int_inds).'; ytar(int_inds).'];

ucomp = chunkerkerneval(chnkr,kerncvel,dens,targs);

utrue = kerns.eval(s,targs)*reshape(strength',[],1);

errplotx= nan(size(xtar));
errploty= nan(size(xtar));

errplotx(int_inds) = abs(ucomp(1:2:end)-utrue(1:2:end));
errploty(int_inds) = abs(ucomp(2:2:end)-utrue(2:2:end));

figure(1)
clf
tiledlayout(1,2,'TileSpacing','compact')
nexttile
h = pcolor(reshape(log10(errplotx),size(xtar))); set(h,'EdgeColor','none');
colormap('jet'); colorbar
clim([-15 -1])
nexttile
h = pcolor(reshape(log10(errploty),size(xtar))); set(h,'EdgeColor','none');
colormap('jet'); colorbar
clim([-15 -1])


%% (Exercise 3) convergence test
% number of chunks to use and number of Legendre nodes on chunks can be
% specified via chunkerfuncuni
%
% - use the above grid of targets for visualizing the error
% - make sure sources are outside (can use the code above)
% - compare true solution (utrue) to computed solution (ucomp) by looping
% over nch = 5*(1:10)
% - compute the error for targets which are inside (use chunkerinterior to
% find targets inside/outside), and store the max of the pointwise l^2 norm
% for each discretization 
% - make a semilogy plot for the error vs nch

% here, we'll check Green's identity using quadarature
% 


% define the combined layer Stokes representation kernel
c = -1;
mu = 1;
kerndvel = kernel('stok','dvel',mu);
kernsvel = kernel('stok','svel',mu);

kerncvel = kerndvel + c*kernsvel;

% convergence
khigh = 12;
bvperrs = zeros(khigh,1);

for k = 1:khigh

  nch = 5*k;
  %
  chnkr = chunkerfuncuni(@(t) starfish(t),nch);
  
  % get boundary data associated to these sources 
  rhs = kerns.eval(s,chnkr)*reshape(strength',[],1);
  
  % get a matrix discretization of the boundary integral equation 
  cmat = chunkermat(chnkr,kerncvel);
  sysmat = cmat - 1/2*eye(2*chnkr.npt) + normonesmat(chnkr);
  
  % solve the system 
  sigma = sysmat\rhs;
  
  % evaluate on a grid
  int_inds = chunkerinterior(chnkr,{xts,yts});
  targs = [];
  targs.r = [xtar(int_inds).'; ytar(int_inds).'];
  
  ucomp = chunkerkerneval(chnkr,kerncvel,sigma,targs);
  
  utrue = kerns.eval(s,targs)*reshape(strength',[],1);
  
  errplotx= nan(size(xtar));
  errploty= nan(size(xtar));
  
  errplotx(int_inds) = abs(ucomp(1:2:end)-utrue(1:2:end));
  errploty(int_inds) = abs(ucomp(2:2:end)-utrue(2:2:end));

  bvperrs(k) = max(sqrt(errplotx(int_inds).^2+errploty(int_inds).^2));
  
  figure(1)
  clf
  tiledlayout(1,2,'TileSpacing','compact')
  nexttile
  h = pcolor(reshape(log10(errplotx),size(xtar))); set(h,'EdgeColor','none');
  colorbar
  clim([-15 -1])
  nexttile
  h = pcolor(reshape(log10(errploty),size(xtar))); set(h,'EdgeColor','none');
  colorbar
  clim([-15 -1])
end

figure(2)
clf
semilogy(5*(1:khigh),bvperrs,'b-o')

%% (Section 4) Stokes flow in a pipe with chunkgraphs
%
% using recursively compressed inverse preconditioning (RCIP) and combined
% layer Stokes representation
%
% we describe problems with corners and multiple junctions using a
% "chunkgraph" structure.  
% 
% NOTE: RCIP in chunkie requires that the diagonal part is the identity. We
% scale below accordingly. 

nedge0 = 4;
xscale = 2*pi;
verts0 = exp(1i*2*pi*(0:3)/4)*exp(1i*pi/4); verts0 = [real(verts0);imag(verts0)];
verts0(1,:) = xscale*verts0(1,:);
vertsx0 = verts0(1,1); vertsy0 = verts0(2,1);
edgesendverts0 = [1:4; 2:4, 1];

fchnks0 = cell(nedge0,1); 
cparams = []; cparams.ta = 0; cparams.tb = 1; cparams.maxchunklen = 0.25;
amp = 0.04; frq = 6*2;
for j = [1 3]
  fchnks0{j} = @(t) sinearc(t,amp,frq);
end

% obstacle in a pipe 
ncurve = 2;
a = -0.5;
b = 0.5;
vertsin = [a b; 0 0];
edgesendvertsin = [1 2; 2, 1];
fchnksin = cell(ncurve,1);
for icurve = 1:ncurve
  fchnksin{icurve} = @(t) [cos(t(:).'); -sin(t(:).')];
end

% outer 1st, inner 2nd
verts = [verts0 vertsin]; % concatenate verts
edgesendverts = [edgesendverts0 size(verts0,2)+edgesendvertsin]; % concatenate edge node info
fchnks = [fchnks0; fchnksin]; % concatenate curves for edges

cg = chunkgraph(verts,edgesendverts,fchnks,cparams);
figure(1),
clf
plot_regions(cg)
hold on
quiver(cg)

%%
% boundary condition specifies the velocity. two values per node
rhs0 =  cg.n ... % direction
      .*[(vertsy0^2-(cg.r(2,:,:)).^2) ... % parabolic profile
      .*(abs(cg.r(2,:,:))<vertsy0) ... % only within vertsy0 limit
      .*(abs(cg.r(1,:,:))>=(vertsx0-1e-10)) ... % only outside vertsx0
      .*sign(cg.r(1,:,:))]; % inlet outlet sign
rhs = rhs0(:); 

figure(1)
clf
plot(cg,'k-')
hold on 
quiver(cg.r(1,:),cg.r(2,:),rhs0(1,:),rhs0(2,:))

%%

% define the combined layer Stokes representation kernel
c = -1;
mu = 1;
kerndvel = -2*kernel('stok','dvel',mu);
kernsvel = -2*kernel('stok','svel',mu);

kerncvel = kerndvel + c*kernsvel;

% get a matrix discretization of the boundary integral equation 
cmat = chunkermat(cg,kerncvel);
sysmat = cmat + eye(2*cg.npt) + normonesmat(merge(cg.echnks));

% solve the system 
sigma = gmres(sysmat,rhs,[],1e-10,100);

% generate some targets...
xtar = linspace(-xscale-1,xscale+1,200);
ytar = linspace(-2,2,200);
[xx,yy] = meshgrid(xtar,ytar);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cg,{xtar,ytar});
uu = nan([2,size(xx)]);

opts = [];
uu(:,in) = reshape(chunkerkerneval(cg,kerncvel,sigma,targs(:,in),opts),2,nnz(in));
u = reshape(uu(1,:,:),size(xx));
v = reshape(uu(2,:,:),size(xx));

figure(2), clf
plot(cg,'k'); hold on
quiver(cg.r(1,:),cg.r(2,:),rhs0(1,:),rhs0(2,:))
quiver(xx,yy,u,v)
%streamslice(xx,yy,u,v)

%% (Exercise 4) 
%
% - add more obstacles in the pipe by correctly specifying verts,
% edgesendverts, and fchnks
% - set the same inlet and outlet flow as the Dirichlet boundary condition
% - solve the new problem and visualize the change of the flow
% 

% outer 1st, inner 2nd
verts = [ verts0 ...                             % outer
          vertsin ...                            % obstacle 1
         [vertsin(1,:)+1; vertsin(2,:)+1/3] ... % obstacle 2
         [vertsin(1,:)-1; vertsin(2,:)-1/3]];   % obstacle 3, concatenate verts
edgesendverts = [edgesendverts0 ...
                 size(verts0,2)+edgesendvertsin ...
                 size(verts0,2)+ncurve+edgesendvertsin ...
                 size(verts0,2)+2*ncurve+edgesendvertsin]; % concatenate edge node info
fchnks = [fchnks0; ...
          fchnksin; ...
          fchnksin; ...
          fchnksin]; % concatenate curves for edges, chunkgraph will snap edges to vertices by translating and rotating

cg = chunkgraph(verts,edgesendverts,fchnks,cparams);

% boundary condition specifies the velocity. two values per node
rhs0 =  cg.n ... % direction
      .*[(vertsy0^2-(cg.r(2,:,:)).^2) ... % parabolic profile
      .*(abs(cg.r(2,:,:))<vertsy0) ... % only within vertsy0 limit
      .*(abs(cg.r(1,:,:))>=(vertsx0-1e-10)) ... % only outside vertsx0
      .*sign(cg.r(1,:,:))]; % inlet outlet sign
rhs = rhs0(:); % interleave

% define the combined layer Stokes representation kernel
c = -1;
mu = 1;
kerndvel = -2*kernel('stok','dvel',mu);
kernsvel = -2*kernel('stok','svel',mu);

kerncvel = kerndvel + c*kernsvel;

% get a matrix discretization of the boundary integral equation 
cmat = chunkermat(cg,kerncvel);
sysmat = cmat + eye(2*cg.npt) + normonesmat(merge(cg.echnks));

% solve the system 
sigma = gmres(sysmat,rhs,[],1e-10,100);

% generate some targets...
xtar = linspace(-xscale-1,xscale+1,200);
ytar = linspace(-2,2,200);
[xx,yy] = meshgrid(xtar,ytar);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cg,{xtar,ytar});
uu = nan([2,size(xx)]);

opts = [];
uu(:,in) = reshape(chunkerkerneval(cg,kerndvel,sigma,targs(:,in),opts),2,nnz(in));
uu(:,in) = uu(:,in) + c*reshape(chunkerkerneval(cg,kernsvel,sigma,targs(:,in),opts),2,nnz(in));
u = reshape(uu(1,:,:),size(xx));
v = reshape(uu(2,:,:),size(xx));

figure(2), clf
plot(cg,'k'); hold on
quiver(cg.r(1,:),cg.r(2,:),rhs0(1,:),rhs0(2,:))
quiver(xx,yy,u,v)
% streamslice(xx,yy,u,v)
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  define functions below here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = stk_singlelayer_eval(s,t,mu)
%STK_SINGLELAYER_EVAL get the interaction matrix for a Stokes single layer
%kernel 
%

% find dimensions 
ns = size(s.r(:,:),2);
nt = size(t.r(:,:),2);

xs = s.r(1,:); ys = s.r(2,:);
xt = t.r(1,:).'; yt = t.r(2,:).';

rx = xt-xs;
ry = yt-ys;

r2 = rx.^2 + ry.^2;

r = sqrt(r2);

log1r = log(1./r);
matxx = 1/(4*pi*mu) * (rx.^2 ./ r2 + log1r);
matyy = 1/(4*pi*mu) * (ry.^2 ./ r2 + log1r);
matxy = 1/(4*pi*mu) * rx.*ry ./ r2;

% Interleave
mat = zeros(2*nt, 2*ns);
mat(1:2:end, 1:2:end) = matxx;
mat(1:2:end, 2:2:end) = matxy;
mat(2:2:end, 1:2:end) = matxy;
mat(2:2:end, 2:2:end) = matyy;

end

function mat = stk_singlelayer_pres_eval(s,t,mu)
%STK_SINGLELAYER_PRES_EVAL get the interaction matrix for a Stokes single
%layer kernel 
%

% find dimensions 
ns = size(s.r(:,:),2);
nt = size(t.r(:,:),2);

xs = s.r(1,:); ys = s.r(2,:);
xt = t.r(1,:).'; yt = t.r(2,:).';

rx = xt-xs;
ry = yt-ys;

r2 = rx.^2 + ry.^2;

matx = 1/(2*pi) * rx ./ r2;
maty = 1/(2*pi) * ry ./ r2;

mat = zeros(nt, 2*ns);
mat(:,1:2:end) = matx;
mat(:,2:2:end) = maty;

end

function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
