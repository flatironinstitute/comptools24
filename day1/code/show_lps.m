% show fund sol and layer pots. Barnett 6/1/24
% needs: arrow.m pkg
clear
figure;
dx = 0.01;
gx1 = 0:dx:2; gx2 = 0:dx:2; [xx1 xx2] = ndgrid(gx1,gx2); y = [1.0,0.9];
rr = sqrt((xx1-y(1)).^2+(xx2-y(2)).^2);
[C,h]=contourf(gx1,gx2,-log(rr)', -0.5:0.5:3.5);
set(h,'linewidth',0.1)
axis equal tight;
fs = 12; xw = 100; yw = 140;
xticks([]); yticks([]);
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
text(y(1)+0.1,y(2)+0.1,'$\mathbf{y}$','interpreter','latex','fontsize',fs);
title('$\Phi(\mathbf{x},\mathbf{y})$','interpreter','latex','fontsize',fs);
%set(gcf,'position',[500 1000 xw yw]); exportgraphics(gcf, '../fs_lap.pdf');

figure
ny = [1 3]; ny = ny/norm(ny);
u = ((xx1-y(1))*ny(1)+(xx2-y(2))*ny(2))./rr.^2;
[C,h]=contourf(gx1,gx2, u', 5*(-1:0.1:1));
set(h,'linewidth',0.1)
arrow(y,y+ny);
axis equal tight;
xticks([]); yticks([]);
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
text(y(1)+ny(1)/2+0.2,y(2)+ny(2)/2,'$n_\mathbf{y}$','interpreter','latex','fontsize',fs);
title('$\partial \Phi(\mathbf{x},\mathbf{y}) / \partial n_\mathbf{y}$','interpreter','latex','fontsize',fs);
%set(gcf,'position',[500 1000 xw yw]); exportgraphics(gcf, '../dipole_lap.pdf','contenttype','vector');

figure;
n=1e3; tt = linspace(0,1,n);
v = [1.5 -0.8]
  z1 = 0.25+v(1)*tt; z2 = 1.5+v(2)*tt;    % param a line
  ny = [-v(2) v(1)]; ny = ny/norm(ny);
slp = 0*xx1; dlp = slp;
for i=1:n
  rr = sqrt((xx1-z1(i)).^2+(xx2-z2(i)).^2);
slp = slp - log(rr) / n;
dlp = dlp + ((xx1-z1(i))*ny(1)+(xx2-z2(i))*ny(2))./rr.^2 / n;
end
[C,h]=contourf(gx1,gx2, slp', -.3 + 0.8*(0:0.1:1));
set(h,'linewidth',0.1)
hold on; plot(z1,z2,'k-','linewidth',2)
axis equal tight;
xticks([]); yticks([]);
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
text(z1(end),z2(end)+0.2,'$\Gamma$','interpreter','latex','fontsize',fs);
set(gcf,'position',[500 1000 xw yw]); exportgraphics(gcf, '../slp_lap.pdf','contenttype','vector');

figure;
[C,h]=contourf(gx1,gx2, dlp', 2*(-1:0.1:1));
set(h,'linewidth',0.1)
hold on; plot(z1,z2,'k-','linewidth',2)
arrow([z1(n/2) z2(n/2)], [z1(n/2)+ny(1), z2(n/2)+ny(2)]);
axis equal tight;
xticks([]); yticks([]);
xlabel('$x_1$','interpreter','latex','fontsize',fs);
ylabel('$x_2$','interpreter','latex','fontsize',fs);
text(z1(end),z2(end)+0.2,'$\Gamma$','interpreter','latex','fontsize',fs);
text(z1(n/2)+ny(1)/2+0.2,z2(n/2)+ny(2)/2,'$n_\mathbf{y}$','interpreter','latex','fontsize',fs);
set(gcf,'position',[500 1000 xw yw]); exportgraphics(gcf, '../dlp_lap.pdf','contenttype','vector');


