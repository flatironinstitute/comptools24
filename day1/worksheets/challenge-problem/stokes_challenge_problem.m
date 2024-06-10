%STOKES_CHALLENGE_PROBLEM
%
% NOTE: This file describes a competition where time is a factor. If you
% plan on competing and the competition is not under way, please stop
% reading now.
%
% Requirements:
%   - this challenge uses the chunkIE package to compute the solution of
%   the BVPs in MATLAB. You can install this package according to the
%   documentation here: chunkie.readthedocs.io. This challenge was designed
%   for v1.0.0 of the chunkIE package. 
%   
% Explanation:
%
% This problem asks you to place 6 obstacles in a pipe in a way that
% minimizes the pressure difference between the two ends for the given flow
% rate. 
%   - you can select the centers and a rotation angle for each of the
%   obstacles (which are starfish shaped)
%   - the obstacles cannot be too close to each other or the pipe walls (no
%   closer than about 0.25 units of separation (about a quarter of the
%   starfish radius)
%   - you must place at least 6 such obstacles (removing them would indeed
%   lower the pressure, so you must put them in)
%   - you can check that the locations meet the requirements by running
%   check_my_geometry (see below for example)
%   - the function get_pressure_difference will set up the boundary value
%   problem and solve and return the pressure difference. It will also
%   return an estimated error for the velocity and pressure difference
%   based on a synthetic test. 
%
% To win the competition, your centers and rotation angles must satisfy the
% geometric constraints. The pressure difference as computed by the 
% original routine get_pressure_difference is what will count (as opposed
% to the value from any modified/sped up solver you might create to find a
% good optimum).
%
% Tips:
% - it is difficult to optimize the center locations automatically because 
% of the constraints. You may want to use some intuition to place these. 
% - the function get_pressure_difference is not optimized for repeated
% calls (indeed it's not exactly optimized for speed). You are encouraged
% to copy this function and make your own version that might be faster
% (perhaps by saving certain re-used quantities, etc.) for testing
% purposes. Note that your answer will be evaluated on the original version
% of the routine.
%

%% pipe set up 

cgpipe = getpipe();
figure(1)
clf
plot_regions(cgpipe)

%% set up and check obstacles (must be at least 6)

% example set up (not a very good one)
ctrs = [-8 -5 5 8 0 0; 0 0 0 0 6 -6]; 
nctr = size(ctrs,2);
rotat = zeros(nctr,1); 

%

[not_in,too_close_pipe,too_close_eachother,cgpipe,chnkrstars,cgpress] = ...
    check_my_geometry(ctrs,rotat);

not_in
too_close_pipe
too_close_eachother

figure(1)
clf
plot(cgpipe,'b-') % pipe
hold on
plot(chnkrstars,'r-') % your obstacles
plot(cgpress,'k--') % where pressure difference is calculated

%%

% last input true tells it to plot 
[press_diff,erru_a,errpress_a] = get_pressure_difference(ctrs,rotat,true)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% put your work here, including any new functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% do not change below here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [press_diff,erru_a,errpress_a] = get_pressure_difference(ctrs,rotat,ifviz)
%GET_PRESSURE_DIFFERENCE
% this sets up the boundary value problem and solves and returns the 
% pressure difference between the front and back pipe ends
%
% Input: 
% ctrs - 2 x nc array of centers for starfish
% rotat - nc x 1 array of angles in radians to rotate starfish
% 
% Optional input:
% ifviz - default(false). if set to true, this will plot some streamlines
%                     for the solution (slow).
%
% Output:
% press_diff - difference of averaged pressure near inlet and outlet of
%                  pipe (quantity to minimize)
% erru_a - relative error in velocity at same points where pressure 
%                  measured for a synthetic test 
% errpress_a - relative error in pressure at same points where pressure 
%                  measured for a synthetic test 
%

if nargin < 3 || isempty(ifviz)
    ifviz = false;
end

assert(size(ctrs,1)==2,'ctrs should be 2 x n array');
nc= size(ctrs(:,:),2);
assert(nc > 5,'must place at least 6 obstacles');

[not_in,too_close_wall,too_close_eachother,cgpipe,chnkrstars,cgpress] = ...
    check_my_geometry(ctrs,rotat);

assert(all(~not_in),'obstacles not in pipe')
assert(all(~too_close_wall),'obstacles too close to wall')
assert(all(~too_close_eachother),'obstacles too close to each other')

% set up BVP

mu = 1; % viscosity
c1 = -2; c2 = 2; % appropriate scalings of single and double for RCIP
kernd = c1*kernel('s','d',mu);
kerns = c2*kernel('s','s',mu);
kernc = kernd + kerns; 

chnkrall = merge([cgpipe.echnks,chnkrstars]); % helps with some tasks to have merged chunker

% pipe flow rhs

vert0 = cgpipe.verts(:,1);
vertsy0 = abs(vert0(2)); vertsx0 = abs(vert0(1));
rhsout = cgpipe.n ... % direction
        .*((vertsy0^2-(cgpipe.r(2,:,:)).^2) ... % parabolic profile
        .*(abs(cgpipe.r(2,:,:))<vertsy0) ... % only within vertsy0 limit
        .*(abs(cgpipe.r(1,:,:))>=(vertsx0-1e-10)) ... % only outside vertsx0
        .*sign(cgpipe.r(1,:,:))); % inlet outlet sign

rhsall1 = zeros(2*chnkrall.npt,1);

rhsall1(1:2*cgpipe.npt) = rhsout(:);

% analytic test rhs

ns = nc + 2;
sources = zeros(2,ns); 
sources(:,1:nc) = ctrs;
sources(:,nc+1) = 0.1*randn(2,1);
sources(:,nc+2) = [-7;7];
src = []; src.r = sources;

strengths = randn(2,ns);

rhsall2 = kerns.eval(src,chnkrall)*strengths(:);

% solve BVPs

sys1 = chunkermat(cgpipe,kernc); % pipe to itself
sys2 = chunkermat(chnkrstars,kernc); % stars to themselves/each other

wts2 = chnkrstars.wts; wts2 = [wts2(:).'; wts2(:).'];
sys12 = kernc.eval(chnkrstars,cgpipe)*diag(wts2(:)); % stars to pipe (constraints make smooth rule accurate)

wts2 = cgpipe.wts; wts2 = [wts2(:).'; wts2(:).'];
sys21 = kernc.eval(cgpipe,chnkrstars)*diag(wts2(:)); % pipe to stars (constraints make smooth rule accurate)

sysall = eye(2*chnkrall.npt) + [sys1 sys12; sys21 sys2] + normonesmat(chnkrall);

dens1 = gmres(sysall,rhsall1,[],1e-10,150);
dens2 = gmres(sysall,rhsall2,[],1e-10,150);

% this kernel set up avoids bugs that users may not have gotten the patch
% for (patch was included in v1.0.0)
kernspres = c2*kernel('s','spres',mu); kernspres.sing = 'pv';
kerndpres = c1*kernel('s','dpres',mu); kerndpres.sing = 'hs';

% analytic error check

utrue = kerns.eval(src,cgpress)*strengths(:);
utest = chunkerkerneval(chnkrall,kernc,dens2,cgpress);

erru_a = norm(utrue-utest)/norm(utrue);

ptrue = kernspres.eval(src,cgpress)*strengths(:);
ptest = chunkerkerneval(chnkrall,kernspres,dens2,cgpress) +  ...
        chunkerkerneval(chnkrall,kerndpres,dens2,cgpress);

avediff = mean(ptrue-ptest); % pressure only determined up to a constant.
errpress_a = norm(ptrue-ptest-avediff)/norm(ptrue);

% pressure for flow problem 

pflow = chunkerkerneval(chnkrall,kernspres,dens1,cgpress) +  ...
        chunkerkerneval(chnkrall,kerndpres,dens1,cgpress);

npt1 = cgpress.echnks(1).npt; % how many points on inflow/outflow lines
npt2 = cgpress.echnks(2).npt; % how many points on inflow/outflow lines
wts1 = cgpress.echnks(1).wts;
wts2 = cgpress.echnks(2).wts;
pflow1 = pflow(1:npt1); pflow2 = pflow((npt1+1):(npt1+npt2));
press_diff = sum(pflow1(:).*wts1(:)) - sum(pflow2(:).*wts2(:));

if ifviz
    % generate some targets...
    xscale = max(chnkrall.r(1,:)); yscale = max(chnkrall.r(2,:));
    x1 = linspace(-xscale+0.1,xscale-0.1,200);
    y1 = linspace(-yscale+0.1,yscale-0.1,200);
    [xx,yy] = meshgrid(x1,y1);
    targs = [xx(:).'; yy(:).'];
    in = chunkerinterior(chnkrall,{x1,y1});
    uu = nan([2,size(xx)]);
    
    uu(:,in) = reshape(chunkerkerneval(chnkrall,kernc,dens1,targs(:,in)),2,nnz(in));
    u = reshape(uu(1,:,:),size(xx));
    v = reshape(uu(2,:,:),size(xx));
    
    figure(2), clf
    plot(chnkrall,'k'); hold on
    quiver(chnkrall.r(1,:).',chnkrall.r(2,:).',rhsall1(1:2:end),rhsall1(2:2:end))
    quiver(xx,yy,u,v)
    streamslice(xx,yy,u,v)
end
   
end

function [not_in,too_close_pipe,too_close_eachother,cgpipe,chnkrstars,cgpress] = ...
    check_my_geometry(ctrs,rotat)
%CHECK_MY_GEOMETRY this function checks that the geometry specified by
%centers and angles is not too close to walls of pipe or each other 
%
% Input: 
% ctrs - 2 x nc array of centers for starfish
% rotat - nc x 1 array of angles in radians to rotate starfish
%
% Output:
% not_in - nc array, not_in(j) true means that jth center not in pipe
% too_close_pipe - nc array, too_close_pipe(j) true means jth starfish too
%    close to pipe wall (must be a quarter radius away, radius approx.
%    0.975
% too_close_eachother - nc array, too_close_eachother(j) true means jth 
%    starfish too close to another starfish
% cgpipe - the pipe geometry (as a chunkgraph)
% chnkrstars - the starfish geometries (as merged chunker)
%

    assert(size(ctrs,1)==2,'ctrs should be 2 x n array');
    nc= size(ctrs(:,:),2);
    assert(numel(rotat)==nc,'number of rotations incompatible with number of centers');

    if nargin < 3
        ifplot = false;
    end

    [cgpipe,cgpress,~,~,~,scal] = getpipe();
    
    not_in = chunkgraphinregion(cgpipe,ctrs).' ~= 2;

    radstar = 1.3*scal;

    xctr = ctrs(1,:); yctr = ctrs(2,:);
    xtot = cgpipe.r(1,:);
    ytot = cgpipe.r(2,:);

    dist_all = sqrt((xtot(:)-xctr).^2 + (ytot(:)-yctr).^2);

    too_close_pipe = min(dist_all,[],1) < radstar*1.25;

    dist_all = sqrt((xctr(:)-xctr).^2 + (yctr(:)-yctr).^2);
    dist_all(1:(nc+1):end) = Inf;

    too_close_eachother = min(dist_all,[],1) < radstar*2.25;

    chnkr0 = chunkerfunc(@(t) starfish(t,3));
    chnkrstars = [];
    for j = 1:size(ctrs,2)
        ctr0 = ctrs(:,j);
        chnkr1 = reverse(chnkr0.move([0;0],ctr0,rotat(j),scal));
        chnkrstars = [chnkrstars,chnkr1];
    end

    chnkrstars = merge(chnkrstars);

end


function [cgpipe,cgpresscheck,R,W,H,scal] = getpipe()
%GETPIPE returns the pipe geometry as a chunkgraph, the surface where the
%pressure is to be measured as a chunkgraph, and parameters of pipe and
%fixes a scaling for the starfish.
%

    % pipe parameters 
    R = 6; W = 4; H = 4;

    % scaling to be used for starfish (default radius is 1.3 so rad will be
    % scal*1.3)
    scal = 0.75;

    verts = [ -(R+W) -R R R+W R+W R -R -(R+W) R/2 -R/2; ...
                H/2  H/2 H/2 H/2 -H/2 -H/2 -H/2 -H/2 0 0 ];
    nv = size(verts,2);
    edgends = [1:8 9 10; 8, 1:7 10 9];

    nedge = size(edgends,2);
    fchnks = cell(nedge,1);
    cparams = cell(nedge,1);
    for j = 1:nedge
        cparams{j}.maxchunklen = scal*1.3;
    end
    
    fchnks{3} = @(t) chnk.curves.bymode(t,3);
    cparams{3}.ta = 0;
    cparams{3}.tb = pi;

    fchnks{7} = @(t) chnk.curves.bymode(t,3);
    cparams{7}.ta = pi;
    cparams{7}.tb = 2*pi;

    fchnks{9} = @(t) [cos(t(:).'); -sin(t(:).')];
    cparams{9}.ta = -pi/8;
    cparams{9}.tb = pi+pi/8;

    fchnks{10} = @(t) [cos(t(:).'); -sin(t(:).')];
    cparams{10}.ta = -pi/8;
    cparams{10}.tb = pi+pi/8;

    cgpipe = chunkgraph(verts,edgends,fchnks,cparams);

    % check pressure just inside the inlet/outlet (save as chunkgraph)

    x = (R+W)-0.2;
    y = H/2-0.2;
    vertspress = [-x, -x, x, x; y -y y -y];
    edgendspress = [1 3; 2 4];
    cgpresscheck = chunkgraph(vertspress,edgendspress);

end