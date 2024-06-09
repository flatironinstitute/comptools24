% Helmholtz exterior Dirichlet BIE solver demo. Agocs 6/8/24

clear
k = 10.0; % wavenumber
%a=0.9; b=1.0;  % 1,0 for unit circle;  0.7, 1 for kite
%Y = @(t) [a*cos(t)+b*cos(2*t); sin(t)];
%Yp = @(t) [-a*sin(t)-2*b*sin(2*t); cos(t)];  % analytic for now
%Ypp = @(t) [-a*cos(t)-4*b*cos(2*t); -sin(t)];   % "
nv = 3;
R = @(t) 1 + 0.3*cos(nv*t);
Rp = @(t) -0.3*nv*sin(nv*t);
Rpp = @(t) -0.3*nv*nv*cos(nv*t);
Y = @(t) [R(t).*cos(t); R(t).*sin(t)];
Yp = @(t) [Rp(t).*cos(t) - R(t).*sin(t); Rp(t).*sin(t) + R(t).*cos(t)];  % analytic for now
Ypp = @(t) [Rpp(t).*cos(t) - 2*Rp(t).*sin(t) - R(t).*cos(t); Rpp(t).*sin(t) + 2*Rp(t).*cos(t) - R(t).*sin(t)];   % "
N = 300;      % (or 300 for _err slide below)
t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
y = Y(t);       % bdry nodes
n = [0 1;-1 0]*Yp(t); speed = sqrt(sum(n.^2,1)); n = n./speed;   % bdry normals
kappa = -sum(Ypp(t) .* n,1)./speed.^2;   % curvatures. check +1 for unit circs

A = zeros(N,N);      % fill A the low-level way...
for i=1:N, for j=1:N
    if j==i
      A(i,j) = -kappa(j)/(2*pi);
    else
      r = y(:,i)-y(:,j);       % r=x-y col vec in R2
      rr = sqrt(sum(r.^2));
      A(i,j) = (1i/2) * k * sum(n(:,j).*r) / rr * besselh(1, k*rr);  % Helm DLP kernel
    end
end, end

A = eye(N) + A*diag(speed.*w);      % note Id doesn't get "speed weights"

% HW is to write the rest of this code:
%f = @(t) 1+0*t;   % test const case, very useful
%f = @(t) cos(5*t+1);   % data vs t, not so common
%ui = @(x) ([1 0]*x) .* ([0 1]*x-0.3);          % test u(x) = x_1(x_2-0.3), not symmetric!
phii = 0;
ui = @(x) exp(1i*k*(cos(phii)*x(1,:) +  sin(phii)*x(2,:))); 
f = @(t) -ui(Y(t));                            % read off its Dirichlet data
rhs = 2*f(t).'; 
sigma = A\rhs;     % solve
rhs(1:10)
sigma(1:10)
                   %sigma = -ones(N,1);  % check DLP[-1] gives +1
figure; dx = 0.01;     % eval solution 
colormap jet;
gx1 = (min(y(1,:))-2):dx:(max(y(1,:))+2); gx2 = (min(y(2,:))-2):dx:(max(y(2,:))+2);  % box it
[xx1 xx2] = ndgrid(gx1,gx2);
u = 0*xx1;
uis = 0*xx1;
for i=1:numel(xx1)        % loop over targ pts
  rr = sqrt(sum(([xx1(i);xx2(i)] - y).^2,1));   % dists from N nodes
  Dker = 1i*k/4*sum(n.*([xx1(i);xx2(i)] - y),1) ./ rr .*besselh(1, k*rr);  % D kernel from nodes
  u(i) = sum(Dker.*speed.*w.*sigma.') + ui([xx1(i); xx2(i)]);     % PTR for D.sigma + incident field

end
mask = inpolygon(xx1,xx2,y(1,:),y(2,:));
imagesc(gx1,gx2, real(u.'), 'alphadata',~mask'); axis xy equal tight;

% R2022b would be needed for alpha in contourf :(
caxis([-1 1.4]); 
caxis_limits = caxis;
max_limit = max(abs(caxis_limits));
caxis([-max_limit max_limit]);
colorbar;

hold on; plot(y(1,:),y(2,:),'k.');   % nodes and l-scaled normals...
%l=0.1; for j=1:N, plot(y(1,j)+[0 l*n(1,j)],y(2,j)+[0 l*n(2,j)],'b-'); end
xticks([]); yticks([]); fs = 12;
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
xtest = [-1.3;2.1];        % one test target far from bdry (far exterior)
plot(xtest(1),xtest(2),'k.','markersize',10);
text(xtest(1)+0.1,xtest(2),'$\mathbf{x}_{test}$','interpreter','latex','fontsize',fs);
text(1.1,0.8,sprintf('N=%d',N));
set(gcf,'position',[500 1000 300 300]); exportgraphics(gcf, '../helmextdir.pdf','contenttype','vector');

%
N = 500;      
t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
y = Y(t);       % bdry nodes
n = [0 1;-1 0]*Yp(t); speed = sqrt(sum(n.^2,1)); n = n./speed;   % bdry normals
kappa = -sum(Ypp(t) .* n,1)./speed.^2;   % curvatures. check +1 for unit circs

A = zeros(N,N);      % fill A the low-level way...
for i=1:N, for j=1:N
    if j==i
      A(i,j) = -kappa(j)/(2*pi);
    else
      r = y(:,i)-y(:,j);       % r=x-y col vec in R2
      rr = sqrt(sum(r.^2));
      A(i,j) = (1i/2) * k * sum(n(:,j).*r) / rr * besselh(1, k*rr);  % Helm DLP kernel
    end
end, end

A = eye(N) + A*diag(speed.*w);      % note Id doesn't get "speed weights"

phii = 0;
ui = @(x) exp(1i*k*(cos(phii)*x(1,:) +  sin(phii)*x(2,:))); 
f = @(t) -ui(Y(t));                            % read off its Dirichlet data
rhs = 2*f(t).'; 
sigma = A\rhs;     % solve
rhs(1:10)
sigma(1:10)
                   %sigma = -ones(N,1);  % check DLP[-1] gives +1
uex = 0*xx1;
uis = 0*xx1;
for i=1:numel(xx1)        % loop over targ pts
  rr = sqrt(sum(([xx1(i);xx2(i)] - y).^2,1));   % dists from N nodes
  Dker = 1i*k/4*sum(n.*([xx1(i);xx2(i)] - y),1) ./ rr .*besselh(1, k*rr);  % D kernel from nodes
  uex(i) = sum(Dker.*speed.*w.*sigma.') + ui([xx1(i); xx2(i)]);     % PTR for D.sigma + incident field
end
% Reference calculation for convergence plot
rtest = sum((xtest - y).^2,1);        % eval at xtest: squared dists from N nodes
Dker = 1i*k/4*sum(n.*(xtest - y),1) ./ rtest .* besselh(1, k*rtest);    % D kernel from nodes
uref = sum(Dker.*speed.*w.*sigma.');         % PTR for D.sigma


figure;
colormap jet;
imagesc(gx1,gx2, log10(abs(u-uex)).', 'alphadata',~mask'); axis xy equal tight;
caxis([-6 0]); colorbar;
hold on; plot(y(1,:),y(2,:),'k.');   % nodes and l-scaled normals...
xticks([]); yticks([]);
text(1.1,0.8,sprintf('N=%d',300));
title('$\log_{10} |u-u_{ref}|$: naive PTR quadr. eval.','interpreter','latex');
set(gcf,'position',[500 1000 300 300]); exportgraphics(gcf, '../helmextdir_err.pdf','contenttype','vector');



%% convergence study...
Ns = 30:30:300; errs = nan*Ns;
for i=1:numel(Ns), N=Ns(i);
  t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
  y = Y(t);       % bdry nodes
  n = [0 1;-1 0]*Yp(t); speed = sqrt(sum(n.^2,1)); n = n./speed;   % bdry normals
  kappa = -sum(Ypp(t) .* n,1)./speed.^2;   % curvatures. check +1 for unit circs
  r1 = y(1,:)'-y(1,:); r2 = y(2,:)'-y(2,:);      % matrix of r=x-y (vec cmpnts 1 & 2)
  rrr = sqrt(r1.^2 + r2.^2);
  A = (1i/2)*k*(n(1,:).*r1 + n(2,:).*r2) ./ rrr .* besselh(1, k*rrr);   % off-diag (-1/pi) n.r/r^2
  A(diagind(A)) = -kappa/(2*pi);                  % overwrite diag elements
  A = eye(N) + A*diag(speed.*w);      % note Id doesn't get "speed weights"
  rhs = 2*f(t).';
  sigma = A\rhs;     % solve
  rr = sum((xtest - y).^2,1);        % eval at xtest: squared dists from N nodes
  Dker = 1i*k/4*sum(n.*(xtest - y),1) ./ rr .* besselh(1, k*rr);    % D kernel from nodes
  utest = sum(Dker.*speed.*w.*sigma.');         % PTR for D.sigma
  errs(i) = (utest-uref);
  fprintf('N=%d\tu err=%.3g\n',N,errs(i))
end
figure; loglog(Ns,abs(errs),'+-'); xlabel('N');
%ylabel('err in $u(\mathbf{x}_{test})$','interpreter','latex');
legend('$|u(\mathbf{x}_{test})-u_{ref}(\mathbf{x}_{test})|$', 'interpreter','latex')
%title('pointwise $N^{-3}$ conv.:','interpreter','latex');
axis([0 max(Ns) 1e-5 10]);
xticks([30 100 300]);
xticklabels({"30", "100", "300"})
set(gcf,'position',[500 1000 200 200]); exportgraphics(gcf, '../helmextdir_conv.pdf','contenttype','vector');
