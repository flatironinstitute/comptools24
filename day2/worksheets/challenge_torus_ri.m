%% Challenge problem - torus
%  
%  In this worksheet, we consider dielectric scattering of a two 
%  tours with major and minor radii 1 and 0.3 respectively, due to 
%  an incident plane wave in the direction [0,0,1] with wavenumber 6\pi.
%
%  The goal of this challenge problem is to find the best refractive
%  index between 1.3 and 1.7 such that energy of the total field due 
%  is maximized in the box  
%  [-0.2, 0.2] \times [-0.2, 0.2] \times [2.1, 3.45]. 
%
%  The index of refraction can be changed by editing line 39.
%
%  You also need to illustrate that your code produces the correct answer
%  by either doing a convergence test, or doing an analytic solution test
%  if you prefer. You may reduce the wavenumber and the refractive index
%  for this part of the exercise if you prefer.
%
%  Utilities for plotting the solution, and computing the energy are
%  provided in this script.

%% Part 1: Verifying order of convergence

run ~/git/fmm3dbie/matlab/startup.m
addpath ~/git/fmm3d/matlab

close('all')
radii = [1,0.3,0];
S = geometries.startorus(radii, 0, [1,1,1], [10,10], 7);
plot(S);

%% Exercise 1: Illustrate accuracy
% Write a test which verifies the solution to the transmission solver

%% Exercise 2: Find the best index of refraction

 zk = 6*pi;

 ri = 1.3; % Find best index of refraction in [1.3, 1.7]
 zks = zk*[1, sqrt(ri)];
 rep_params = ones(4,1);

 dir = [0,0,1];
 [uin, uingrad] = helm3d.planewave(zks(1), dir, S);
 
 rhs = complex(zeros(2,S.npts));
 rhs(1,:) = -uin;
 rhs(2,:) = -(uingrad(1,:).*S.n(1,:) + uingrad(2,:).*S.n(2,:) + ...
     uingrad(3,:).*S.n(3,:));

 eps = 1e-4;
 opts = [];
 opts.eps_gmres = 1E-6;
 t1 = tic; [sig, errs] = solver(S, 'helm', 'trans', rhs, eps, zks, ...
                           rep_params,opts);
 tend = toc(t1);
 fprintf('Time taken in solver=%d\n',tend);

 %% Evaluate field values on a slice

xmin = -0.2;
xmax =  0.2;

ymin = -0.2;
ymax = 0.2;

zmin = 2.1;
zmax = 3.45;

blims2d = [xmin, zmin; xmax, zmax];
xsec = [0,1,0];
xsec = xsec.'/norm(xsec);

[X, Y, Z] = get_2dslice(blims2d, xsec); 

sz = size(X);
xyz_out = [X(:).';Y(:).';Z(:).'];

targ_info = [];
targ_info.r = xyz_out;

t2 = tic;
eps = 1E-4;
pot = eval_fields(S, 'h', 'trans', sig, targ_info, eps, zks, rep_params); 
tend = toc(t2);
fprintf('time taken in eval=%d\n', tend)

pot_in = helm3d.planewave(zks(1), dir, targ_info);
utot = pot_in + pot;
utot_plt = reshape(utot, sz);

%% Some plotting utilities to guide your search

figure;
h = plot_surface_with_slice(S, X, Y, Z, abs(utot_plt));


%% Compute the energy

klege = 20;
blims3d = [xmin, ymin, zmin; xmax, ymax, zmax];
[X, Y, Z] = get_3dvol(klege, blims3d);

targs = [];
targs.r = [X.';Y.';Z.'];
t2 = tic;
eps = 1E-6;
pot = eval_fields(S, 'h', 'trans', sig, targs, eps, zks, rep_params); 
tend = toc(t2);

pot_in = helm3d.planewave(zks(1), dir, targs);
utot = pot_in + pot;

energy = sum(abs(utot).^2.*W);






function [X, Y, Z] = get_2dslice(blims, xsec)
%
% Create a 2d slice in 3d
%
% Input arguments:
%   - blims: [umin, vmin; umax, vmax] defining the limits of the plane
%   - xsec: normal to the plane
%

    umin = blims(1,1);
    umax = blims(2,1);

    vmin = blims(1,2);
    vmax = blims(2,2);

    vnulls = null(xsec.');
    uaxi = vnulls(:,1);
    vaxi = vnulls(:,2);

    us = umin:0.05:umax;
    vs = vmin:0.05:vmax;
    [U,V] = meshgrid(us,vs);

    X = U*uaxi(1) + V*vaxi(1);
    Y = U*uaxi(2) + V*vaxi(2);
    Z = U*uaxi(3) + V*vaxi(3);

end



function [X, Y, Z, W] = get_3dvol(k, blims)

    [x, w, ~] = polytens.lege.pts(k);
    [X, Y, Z] = ndgrid(x,x,x);

    xmin = blims(1,1);
    xmax = blims(2,1);

    ymin = blims(1,2);
    ymax = blims(2,2);
    
    zmin = blims(1,3);
    zmax = blims(2,3);

    [WX, WY, WZ] = ndgrid(w,w,w);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    WX = WX(:);
    WY = WY(:);
    WZ = WZ(:);
    X = (X + 1)/2*(xmax - xmin) + xmin;
    Y = (Y + 1)/2*(ymax - ymin) + ymin;
    Z = (Z + 1)/2*(zmax - zmin) + zmin;
    
    WX = WX*(xmax - xmin)/2;
    WY = WY*(ymax - ymin)/2;
    WZ = WZ*(zmax - zmin)/2;
    W = WX.*WY.*WZ;


end


function h = plot_surface_with_slice(S, X, Y, Z, f)
% Plot surface with volumetric slice of data
% 
% Input arguments:
%   S - surfer object
%   X, Y, Z - meshgrid of slice data
%   f - function values on slice (must be real)
%  
% Output arguments:
%   h - object handle for plot
%

    h = surf(X, Y, Z, f, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.9);
    colorbar;
    lims = clim;
    
    hold on
    
    h = plot(S);
    set(h,'EdgeColor','black');
    set(h,'EdgeAlpha',0.3);
    set(h,'FaceColor','none');
    caxis(lims);
end
