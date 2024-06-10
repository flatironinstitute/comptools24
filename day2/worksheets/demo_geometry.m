% Worksheet 1: Geometry and integrals
%
%  Learning Goals:
%    - Discretizing surfaces defined via analytic parametrizations
%    - Plotting, and manipulating these surfaces,
%        including (moving, rotating, and merging them)
%    - Understanding the data structure used to represent the surface
%    - Computing integrals of smooth functions on surfaces
%    - Estimating error in functions discretized on surface
%    - Oversampling the discretized surface
%    - Inputting high-order geometries via smoothing low-order ones
%      using the surface smoother.
%  
%
%% Section 1: Load simple geometries (Step 3 of the 6 step program)
%
%
clearvars 
clear all
clear classes

%% Run startup script
run ~/git/fmm3dbie/matlab/startup.m %Enter your path to fmm3dbie startup ehre

%%
close('all')
R = 1;
S = geometries.sphere(R, 4);

% Plot the domain
figure
plot(S);

% Compute errors
figure
errps = surf_fun_error(S, S.mean_curv);
errps = log10(errps);
plot(S, errps);

%% 
% scatter plot
close('all')
figure
scatter(S);
axis equal

figure
plot_nodes(S);

%%
% Look inside S 
S

%% Section 2: Computing integrals and error estimatation 
% Compute the surface area, centroid, and moment of inertia
a = sum(S.wts(:));
erra = abs(a - 4*pi*R^2);
fprintf('Error in computing surface area of sphere=%d\n',erra);


%% Exercise 1: 
% Use this to compute the centroid, and moment of interia 


%% Harder integral (Gauss' Law)
% \int_{\Gamma} \nabla_{x} \frac{1}{4 \pi |x-y|} \cdot n(x) dx = \delta_{S}

close('all')
xyz = -[0.3; 0.2; 0.1];
dx = S.r(1,:) - xyz(1);
dy = S.r(2,:) - xyz(2);
dz = S.r(3,:) - xyz(3);
r = sqrt(dx.^2 + dy.^2 + dz.^2);

rdotn = dx.*S.n(1,:) + dy.*S.n(2,:) + dz.*S.n(3,:);
fint = (rdotn./r.^3/4/pi);

% Plot surface with function
figure
plot(S, fint)

errps = surf_fun_error(S, fint);
errps = log10(errps);
figure
plot(S, errps)

err = (fint*S.wts - 1);
fprintf('Error at easy point = %d\n', err);

%% Harder test, using anonymous functions
% Move point closer to boundary -> make a function
close('all')
xyz = -[0.95; 0.01; 0.03]*R;
dx = @(S, xyz) S.r(1,:) - xyz(1);
dy = @(S, xyz) S.r(2,:) - xyz(2);
dz = @(S, xyz) S.r(3,:) - xyz(3);
r = @(S, xyz) sqrt(dx(S, xyz).^2 + dy(S, xyz).^2 + dz(S, xyz).^2);

rdotn = @(S, xyz) dx(S, xyz).*S.n(1,:) + ...
                   dy(S, xyz).*S.n(2,:) + dz(S, xyz).*S.n(3,:);
fint = @(S, xyz) (rdotn(S, xyz)./r(S, xyz).^3/4/pi);

fint1 = fint(S, xyz);

plot(S, fint1);

errps = surf_fun_error(S, fint1);
errps = log10(errps);
figure
plot(S, errps)

err = (fint1*S.wts - 1);
fprintf('Error at difficult point = %d\n', err);

%% Section 3: Oversampling
S2 = oversample(S, 20);

fint2 = fint(S2, xyz);

errps = surf_fun_error(S2, fint2);
errps = log10(errps);
figure
plot(S2, errps)

err = (fint2*S2.wts - 1);
fprintf('Error at difficult point post oversampling= %d\n', err);

%% Ellipsoid geometries
close('all')
abc = [2;3.1; 1.7];
S_ellip = geometries.ellipsoid(abc, [2,2,2], [0,0,0], 6);
plot(S_ellip)


%% Exercise 2/Section 4: Manipulate geometries
% Verify area using mathematica/wikipedia or
% elliptic integrals if you are brave, and the repeat the gauss test

% Rotate

% translate

% scale


%% Section 5: Surface smoother
% Surface smoother (supports many different low order meshes)

opts.nrefine = 2;
opts.rlam = 2.5;
S = multiscale_mesher('../geometries_flat/simple_torus.msh', 6, opts);

%%
close('all')
figure
plot(S{2})

errps = surf_fun_error(S{2}, S{2}.mean_curv);
figure
plot(S{2}, log10(errps));

figure
plot(S{3})

errps = surf_fun_error(S{3}, S{3}.mean_curv);
figure
plot(S{3}, log10(errps));

%% User exercise: Play around with smoothing parameter, and a different
% geometry, plot the errors in refined surfaces, see convergence


