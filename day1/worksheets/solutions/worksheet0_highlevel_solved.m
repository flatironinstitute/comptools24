%WORKSHEET0_HIGHLEVEL_SOLVED
%
% NOTE: these are the solutions of the worksheet! It's recommended to try
% the worksheet without consulting these solutions first.
%
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - this worksheet goes over an integral equation solve using the chunkIE
%   software. 
%   - At a high level, the user creates the geometry, selects an
%   appropriate kernel, builds the system matrix for the BIE, solves the
%   system matrix for the BIE, and then plots the result
%   - Shows a convergence test
%   - Also demonstrates a problem with corners 
%

%% Create a geometry, and plot it

chnkr = chunkerfunc(@(t) starfish(t)); % chunker object from function handle

% always good to check the normals
figure(1)
clf
plot(chnkr)
hold on
quiver(chnkr);

%% Let's try an interior Dirichlet problem for laplace first
%
% with the representation u = D[sigma], the boundary integral equation 
% is -1/2 sigma + D[sigma] = f, where f is the Dirichlet boundary data

kern = kernel('lap','d'); % kernel object 
sysmat = chunkermat(chnkr,kern);

% and add the appropriate multiple of the identity
sysmat = -0.5*eye(chnkr.npt)+sysmat;

% we can manufacture an analytic solution x^2-y^2 

fun = @(x,y) x.^2-y.^2;

x_geo = chnkr.r(1,:).';
y_geo = chnkr.r(2,:).';
rhs = fun(x_geo,y_geo);

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

[fints] = chunkerkerneval(chnkr,kern,dens,targs);

soln = NaN*ones(sz);
soln(int_inds) = fints;

soltrue = NaN*soln;
soltrue(int_inds) = fun(xtar(int_inds),ytar(int_inds));

err = abs(soln - soltrue)/max(abs(soltrue(:)));

figure(1)
clf
h = pcolor(xts,yts,reshape(log10(abs(err)),sz));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'k-')
colorbar
axis equal

%% Convergence

nstep = 8;
soln = zeros(nrefs,1);

src = [1;-1];
tar = [0;0.1];
targs = []; 
targs.r = tar;

for ii=1:nstep

    nch = ii*5;
    chnkr2 = chunkerfuncuni(@(t) starfish(t),nch);
    kern = kernel('lap','d');
    sysmat = chunkermat(chnkr2,kern);

    % and add the appropriate multiple of the identity
    sysmat = -0.5*eye(chnkr2.npt)+sysmat;

    rhs = chnk.lap2d.green(src,chnkr2.r(:,:));

    % and solve
    dens = sysmat\rhs;

    [fints] = chunkerkerneval(chnkr2,kern,dens,targs);
    soln(ii) = fints;
end

soltrue = chnk.lap2d.green(src,tar);

abs(soln-soltrue)


%% (Exercise 1) 
%
% devise your own analytical solution of the laplace equation and re-run
% the steps above
%
% advanced option: solve the corresponding Helmholtz Dirichlet problem by
% selecting a wave number zk and getting the kernel 
%           kernd = kernel('h','d',zk).
% can you devise an analytic solution of (Delta + k^2)u = 0
% note: if zk has large real part, might need more chunks. see
% help chunkerfunc or help chunkerfuncuni
%

% laplace part 

fun = @(x,y) exp(x).*cos(y);

%% 
% helmholtz part

zk = 2.1;

cparams = []; cparams.maxchunklen = 4.0/zk;
chnkr = chunkerfunc(@(t) starfish(t),cparams);

kerndh = kernel('helm','d',zk); % kernel object 
sysmat = chunkermat(chnkr,kerndh);

% and add the appropriate multiple of the identity
sysmat = -0.5*eye(chnkr.npt)+sysmat;

% we can manufacture an analytic solution 
funhelm = @(x,y,zk) cos(zk*x) + sin(zk*y);

x_geo = chnkr.r(1,:).';
y_geo = chnkr.r(2,:).';
rhs = funhelm(x_geo,y_geo,zk);

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

[fints] = chunkerkerneval(chnkr,kerndh,dens,targs);

soln = NaN*ones(sz);
soln(int_inds) = fints;

soltrue = NaN*soln;
soltrue(int_inds) = funhelm(xtar(int_inds),ytar(int_inds),zk);

err = abs(soln - soltrue)/max(abs(soltrue(:)));

figure(2)
clf
h = pcolor(xts,yts,reshape(log10(abs(err)),sz));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'k-');
colorbar
axis equal
  
