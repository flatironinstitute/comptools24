% Nystrom method demos for SKIE on periodic interval. Barnett 5/30/24
clear

if 1 % 1) analytically known test ----------------------------------------------
kfun = @(t,s) exp(3*cos(t-s));   % smooth, convolutional kernel, domain [0,2pi)
ffun = @(t) cos(5*t+1);      % data (RHS) func
sigexfun = @(t) cos(5*t+1) / (1 + 2*pi*besseli(5,3));   % soln known
% Homework question: show this is the exact solution [Hints: diag in Fourier basis, and look up integral form for I_nu(z) Bessel func]

Ns = 10:5:40;            % convergence study
sigtests=0*Ns;
for i=1:numel(Ns), N=Ns(i);
  t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % nodes, weights, row vecs
  A = eye(N) + bsxfun(kfun,t',t)*diag(w);  % fills k(t_i,t_j) w_j for i,j=1..N
  rhs = ffun(t');          % col vec
  sig = A\rhs;             % dense direct solve (pivoted LU)
  sigtests(i) = sig(end);      % value of soln func sigma(0), always last node
  fprintf("%d\t%.12g\n", N, sigtests(i))
end

figure; subplot(1,3,1); imagesc(A);
axis equal tight; title('system matrix A'); colorbar;
subplot(1,3,2); tt = linspace(0,2*pi,1e3);   % plot grid
plot(tt, ffun(tt), 'b-', t, sig, 'k+', tt, sigexfun(tt), 'g-');
xlabel('t'); ylabel('funcs'); legend('RHS func f(t)', 'soln vec \sigma_j', 'exact \sigma(t)');
errs = abs(sigtests-sigtests(end));   % estimated errors ("self convergence")
subplot(1,3,3); semilogy(Ns,errs,'+-'); xlabel('N'); ylabel('err in \sigma(0)');
errsex = abs(sigtests-sigexfun(0));   % actual errors
hold on; semilogy(Ns,errsex,'ro-'); axis tight; legend('self-conv','true err');
set(gcf,'position',[500 1000 1000 250]); exportgraphics(gcf, '../nyst_conv.pdf');
end




if 1 % 2) demo some reduced conv rates -----------------------------------------
  % 2a)
kfun = @(s,t) exp(3*cos(t-s));
%kfun = @(s,t) exp(2*cos(t-s+2*sin(s)));   % smooth general kernel, domain [0,2pi)
ffun = @(t) abs(sin(t));      % data (RHS) func, C not C^1, kink not at origin

Ns = 40:40:1000;            % convergence study
sigtests=0*Ns;
for i=1:numel(Ns), N=Ns(i);
  t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % nodes, weights, row vecs
  A = eye(N) + bsxfun(kfun,t',t)*diag(w);  % fills k(t_i,t_j) w_j for i,j=1..N
  rhs = ffun(t');         % col vec
  sig = A\rhs;             % dense direct solve (pivoted LU)
  sigtests(i) = sig(end);      % value of soln func sigma(0), always last node
  fprintf("%d\t%.12g\n", N, sigtests(i))
end
figure;
errs = abs(sigtests-sigtests(end));   % estimated errors ("self convergence")
subplot(1,2,1); loglog(Ns,errs,'+-', Ns,5*Ns.^(-2),'r--');
axis tight; xlabel('N'); legend('self-conv err in $\sigma(0)$','$1/N^2$','interpreter','latex');
fs = 12;
text(50,1e-5,'$f(t) = |\sin t|$','interpreter','latex','fontsize',fs);
text(50,2e-6,'$k$ smooth','interpreter','latex','fontsize',fs);

% 2b)
kfun = @(s,t) 10*sin(abs(t-s)/2).^3;  % nonsmooth, conv kernel
% Homework: is k continuous? Is the resulting K compact?
ffun = @(t) cos(5*t+1);      % data (RHS) func, smooth again
Ns = 20:20:500;            % convergence study
sigtests=0*Ns;
for i=1:numel(Ns), N=Ns(i);
  t = 2*pi/N*(1:N); w = 2*pi/N*ones(1,N);   % nodes, weights, row vecs
  A = eye(N) + bsxfun(kfun,t',t)*diag(w);  % fills k(t_i,t_j) w_j for i,j=1..N
  rhs = ffun(t');          % col vec
  sig = A\rhs;             % dense direct solve (pivoted LU)
  sigtests(i) = sig(end);      % value of soln func sigma(0), always last node
  fprintf("%d\t%.12g\n", N, sigtests(i))
end
errs = abs(sigtests-sigtests(end));   % estimated errors ("self convergence")
subplot(1,2,2); loglog(Ns,errs,'+-', Ns,5*Ns.^(-4),'r--');
axis tight; xlabel('N'); legend('self-conv err in $\sigma(0)$','$1/N^4$','interpreter','latex');
text(25,3e-9,'$f$ smooth','interpreter','latex','fontsize',fs);
text(25,3e-10,'$k(t,s)=\sin^3 \frac{|t-s|}{2}$','interpreter','latex','fontsize',fs);  % ignore prefac

set(gcf,'position',[500 1000 600 240]); exportgraphics(gcf, '../nyst_nonsm.pdf');
end
