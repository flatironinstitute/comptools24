% Simplest Laplace interior Dirichlet BIE solver demo. Barnett 6/3/24
clear
a=0.7; b=1.0;  % 1,0 for unit circle;  0.7, 1 for kite
y = @(t) [a*cos(t)+b*cos(2*t); sin(t)];
yp = @(t) [-a*sin(t)-2*b*sin(2*t); cos(t)];  % analytic for now
ypp = @(t) [-a*cos(t)-4*b*cos(2*t); -sin(t)];   % "
N = 100;
tj = 2*pi/N*(1:N); wj = 2*pi/N*ones(1,N);   % PTR nodes & weights, row vecs
yj = y(tj);       % bdry nodes
nj = [0 1;-1 0]*yp(tj); speedj = sqrt(sum(nj.^2,1)); nj = nj./speedj;   % bdry normals
kappaj = -sum(ypp(tj) .* nj,1)./speedj.^2;   % curvatures. check +1 for unit circs
A = zeros(N,N);
for i=1:N, for j=1:N
    if j==i
      A(i,j) = kappaj(j)/(2*pi);
    else
      d = yj(:,i)-yj(:,j);       % x-y col vec in R2
      A(i,j) = (-1/pi) * sum(nj(:,j).*d) / sum(d.^2);  % Lap DLP kernel
    end
end, end
% or, a more "idomatic matlab" way to fill A...
%dd1 = yj(1,:)'-yj(1,:); dd2 = yj(2,:)'-yj(2,:);  % 2 cmpnts of x-y matrix
%A = (-1/pi)*(nj(1,:).*dd1 + nj(2,:).*dd2) ./ (dd1.^2+dd2.^2); 
%A(diagind(A)) = kappaj/(2*pi);     % a utility
A = eye(N) + A*diag(speedj.*wj);    % note Id doesn't get "speed weights"
                                    %f = @(t) 1+0*t;
f = @(t) cos(5*t+1);   % data
rhs = -2*f(tj)';
sigma = A\rhs;     % solve
                   %sigma = -ones(N,1);  % check DLP[-1] gives +1
figure; dx = 0.01;     % eval solution 
gx1 = min(yj(1,:)):dx:max(yj(1,:)); gx2 = min(yj(2,:)):dx:max(yj(2,:));
[xx1 xx2] = ndgrid(gx1,gx2);
u = 0*xx1;
for i=1:numel(xx1)        % loop over targ pts
  r2 = sum(([xx1(i);xx2(i)] - yj).^2,1);   % squared dists from N nodes
  Dkerj = sum(nj.*([xx1(i);xx2(i)] - yj),1) ./ (2*pi*r2);  % D kernel from nodes
  u(i) = sum(Dkerj.*speedj.*wj.*sigma');     % PTR for D.sigma
end
mask = inpolygon(xx1,xx2,yj(1,:),yj(2,:));
surf(gx1,gx2, 0*u', u', 'alphadata',mask');
caxis([0 2]); colorbar;
hold on; plot(yj(1,:),yj(2,:),'k.');
l=0.1; for j=1:N, plot(yj(1,j)+[0 l*nj(1,j)],yj(2,j)+[0 l*nj(2,j)],'b-'); end
axis equal tight;
xticks([]); yticks([]);
fs = 12;
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
set(gcf,'position',[500 1000 300 300]); exportgraphics(gcf, '../lapintdir.pdf','contenttype','vector');
