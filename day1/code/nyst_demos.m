% Nystrom method demos for SKIE on periodic interval. Barnett 5/30/24

clear

% analytically known solution:
kfun = @(s,t) exp(3*cos(t-s));   % smooth, convolutional kernel, domain [0,2pi)
ffun = @(t) cos(5*t+1);      % data (RHS) func
sigexfun = @(t) cos(5*t+1) / (1 + 2*pi*besseli(5,3));   % analytic soln

Ns = 10:5:50;            % convergence study
sigtests=0*Ns;
for i=1:numel(Ns), N=Ns(i);
  tj = 2*pi/N*(1:N);       % nodes, row vec
  wj = 2*pi/N*ones(1,N);   % weights, row vec
  A = eye(N) + bsxfun(kfun,tj',tj) .* wj;  % fills K(i,j) w_j for i,j=1..N
                                           % note I doens't get weighted
  rhs = ffun(tj');         % col vec
  sig = A\rhs;             % dense direct solve (pivoted LU)
  sigtests(i) = sig(end);      % value of soln func sigma(0), always last node
  fprintf("%d\t%.12g\n", N, sigtests(i))
end

figure; subplot(1,3,1); imagesc(A);
axis equal tight; title("system matrix A"); colorbar;
subplot(1,3,2); tt = linspace(0,2*pi,1e3);
plot(tt, ffun(tt), 'b-', tj, sig, 'k+', tt, sigexfun(tt), 'g-');
xlabel('t'); ylabel('funcs'); legend('f(t)', '\sigma_j', 'exact \sigma(t)');
errs = abs(sigtests-sigtests(end));   % estimated errors ("self convergence")
subplot(1,3,3); semilogy(Ns,errs,'+-'); xlabel("N"); ylabel("err in \sigma(0)");

