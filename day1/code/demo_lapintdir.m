% Simplest Laplace interior Dirichlet BIE solver demo. Barnett 6/3/24
% Can generate 3 figs for slides.

clear
a=0.9; b=1.0;  % 1,0 for unit circle;  0.7, 1 for kite
Y = @(t) [a*cos(t)+b*cos(2*t); sin(t)];
Yp = @(t) [-a*sin(t)-2*b*sin(2*t); cos(t)];  % analytic for now
Ypp = @(t) [-a*cos(t)-4*b*cos(2*t); -sin(t)];   % "
N = 100;      % (or 300 for _err slide below)
t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
y = Y(t);       % bdry nodes
n = [0 1;-1 0]*Yp(t); speed = sqrt(sum(n.^2,1)); n = n./speed;   % bdry normals
kappa = -sum(Ypp(t) .* n,1)./speed.^2;   % curvatures. check +1 for unit circs

A = zeros(N,N);      % fill A the low-level way...
for i=1:N, for j=1:N
    if j==i
      A(i,j) = kappa(j)/(2*pi);
    else
      r = y(:,i)-y(:,j);       % r=x-y col vec in R2
      A(i,j) = (-1/pi) * sum(n(:,j).*r) / sum(r.^2);  % Lap DLP kernel
    end
end, end
% or, a more "idomatic matlab" way to fill A...
r1 = y(1,:)'-y(1,:); r2 = y(2,:)'-y(2,:);      % matrix of r=x-y (vec cmpnts 1 & 2)
A = (-1/pi)*(n(1,:).*r1 + n(2,:).*r2) ./ (r1.^2+r2.^2);   % off-diag (-1/pi) n.r/r^2
A(diagind(A)) = kappa/(2*pi);                  % overwrite diag elements

A = eye(N) + A*diag(speed.*w);      % note Id doesn't get "speed weights"

% HW is to write the rest of this code:
%f = @(t) 1+0*t;   % test const case, very useful
%f = @(t) cos(5*t+1);   % data vs t, not so common
uex = @(x) ([1 0]*x) .* ([0 1]*x-0.3);          % test u(x) = x_1(x_2-0.3), not symmetric!
f = @(t) uex(Y(t));                            % read off its Dirichlet data
rhs = -2*f(t)';
sigma = A\rhs;     % solve
                   %sigma = -ones(N,1);  % check DLP[-1] gives +1
figure; dx = 0.01;     % eval solution 
gx1 = min(y(1,:)):dx:max(y(1,:)); gx2 = min(y(2,:)):dx:max(y(2,:));  % box it
[xx1 xx2] = ndgrid(gx1,gx2);
u = 0*xx1;
for i=1:numel(xx1)        % loop over targ pts
  rr = sum(([xx1(i);xx2(i)] - y).^2,1);   % squared dists from N nodes
  Dker = sum(n.*([xx1(i);xx2(i)] - y),1) ./ (2*pi*rr);  % D kernel from nodes
  u(i) = sum(Dker.*speed.*w.*sigma');     % PTR for D.sigma
end
mask = inpolygon(xx1,xx2,y(1,:),y(2,:));
imagesc(gx1,gx2, u', 'alphadata',mask'); axis xy equal tight;
% R2022b would be needed for alpha in contourf :(
caxis([-1 1.4]); colorbar;
hold on; plot(y(1,:),y(2,:),'k.');   % nodes and l-scaled normals...
%l=0.1; for j=1:N, plot(y(1,j)+[0 l*n(1,j)],y(2,j)+[0 l*n(2,j)],'b-'); end
xticks([]); yticks([]); fs = 12;
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
xtest = [0.8;0.3];        % one test target far from bdry (deep interior)
plot(xtest(1),xtest(2),'k.','markersize',10);
text(xtest(1)+0.1,xtest(2),'$\mathbf{x}_{test}$','interpreter','latex','fontsize',fs);
text(1.1,0.8,sprintf('N=%d',N));
set(gcf,'position',[500 1000 300 300]); exportgraphics(gcf, '../lapintdir.pdf','contenttype','vector');


if 0   % rerun above at N=300 for this error plot....
figure; uuex = reshape(uex([xx1(:),xx2(:)]'), size(xx1));
imagesc(gx1,gx2, log10(abs(u-uuex))', 'alphadata',mask'); axis xy equal tight;
caxis([-13 0]); colorbar;
hold on; plot(y(1,:),y(2,:),'k.');   % nodes and l-scaled normals...
xticks([]); yticks([]);
text(1.1,0.8,sprintf('N=%d',N));
title('$\log_{10} |u-u_{ex}|$: naive PTR quadr. eval.','interpreter','latex');
set(gcf,'position',[500 1000 300 300]); exportgraphics(gcf, '../lapintdir_err.pdf','contenttype','vector');
end



% convergence study...
Ns = 30:30:240; errs = nan*Ns;
for i=1:numel(Ns), N=Ns(i);
  t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
  y = Y(t);       % bdry nodes
  n = [0 1;-1 0]*Yp(t); speed = sqrt(sum(n.^2,1)); n = n./speed;   % bdry normals
  kappa = -sum(Ypp(t) .* n,1)./speed.^2;   % curvatures. check +1 for unit circs
  r1 = y(1,:)'-y(1,:); r2 = y(2,:)'-y(2,:);      % matrix of r=x-y (vec cmpnts 1 & 2)
  A = (-1/pi)*(n(1,:).*r1 + n(2,:).*r2) ./ (r1.^2+r2.^2);   % off-diag (-1/pi) n.r/r^2
  A(diagind(A)) = kappa/(2*pi);                  % overwrite diag elements
  A = eye(N) + A*diag(speed.*w);      % note Id doesn't get "speed weights"
  rhs = -2*f(t)';
  sigma = A\rhs;     % solve
  rr = sum((xtest - y).^2,1);        % eval at xtest: squared dists from N nodes
  Dker = sum(n.*(xtest - y),1) ./ (2*pi*rr);    % D kernel from nodes
  utest = sum(Dker.*speed.*w.*sigma');         % PTR for D.sigma
  errs(i) = utest-uex(xtest);
  fprintf('N=%d\tu err=%.3g\n',N,errs(i))
end
figure; semilogy(Ns,abs(errs),'+-'); xlabel('N');
%ylabel('err in $u(\mathbf{x}_{test})$','interpreter','latex');
legend('$|u(\mathbf{x}_{test})-u_{ex}(\mathbf{x}_{test})|$', 'interpreter','latex')
title('pointwise spectral conv.:','interpreter','latex');
axis([0 max(Ns) 1e-13 1]);
set(gcf,'position',[500 1000 200 200]); exportgraphics(gcf, '../lapintdir_conv.pdf','contenttype','vector');
