clear all; close all; clc

%% Initial Definitions
tspan = 0:1:200;
nx = 64;
ny = nx;
n = nx*ny;
xspan = [-10 10];
Lx = xspan(2)-xspan(1);
delta_x = (xspan(2)-xspan(1))/nx;
yspan = [-10 10];
Ly = yspan(2)-yspan(1);
delta_y = (yspan(2)-yspan(1))/ny;
xgrid = linspace(xspan(1),xspan(2),nx+1);
xgrid = xgrid(1:end-1);
ygrid = linspace(yspan(1),yspan(2),ny+1);
ygrid = ygrid(1:end-1);
[X,Y] = meshgrid(xgrid,ygrid);

% u_0 = exp(-X.^2/3-Y.^2);
% u_0 = (cos(X).^2).*exp(-Y.^2);
% u_0 = reshape(u_0,n,1);
% v_0 = exp(-X.^2-Y.^2/3);
% v_0 = reshape(v_0,n,1);


u_0 = normrnd(1,0.15,n,1);
v_0 = normrnd(1,0.15,n,1);

% u_0 = ones(size(X));
% u_0((abs(X) < Lx/4) & (abs(Y) < Ly/4)) = 0.5;
% u_0 = u_0 + normrnd(0,0.0075,size(u_0));
% u_0 = reshape(u_0,n,1);
% 
% v_0 = zeros(size(X));
% v_0((abs(X) < Lx/4) & (abs(Y) < Ly/4)) = 0.25;
% v_0 = v_0 + normrnd(0,0.0075,size(v_0));
% v_0 = reshape(v_0,n,1);

uv_0 = [u_0; v_0];

ru = 0.00002;
rv = 0.00001;
f = 0.04;
k = 0.06;

%% Build derivative matrices
 
 %build dy

ov = ones(nx,1);
sub_d = spdiags([ov -ov ov -ov], [-(nx-1) -1 1 (nx-1)],nx,nx); %nx*nx matrix appearing along the diagonal of the full n*n
sub_od = zeros(nx); %nx*nx matrix appearing off the diagonal of the full n*n

d_y = zeros(n);
for j = 1:nx
    for k = 1:nx
        if j==k
            d_y((j-1)*nx+1:j*nx,(k-1)*nx+1:k*nx) = sub_d;
        else
            d_y((j-1)*nx+1:j*nx,(k-1)*nx+1:k*nx) = sub_od;
        end
    end
end
d_y = d_y/(2*delta_y);
d_y = sparse(d_y);

%build dx

ov = ones(n,1);
d_x = spdiags([ov -ov ov -ov], [-(n-nx) -nx nx (n-nx)],n,n);
d_x = d_x/(2*delta_x);

% build Laplacian
zv = zeros(n,1);
e0 = zv;
e1 = ov;
e2 = ov;
e4 = zv;

for j = 1:nx
    e2(nx*j) = 0;
    e4(nx*j) = 1;
end
e3(2:n,1) = e2(1:n-1,1);
e3(1,1) = e2(n,1);
e5(2:n,1) = e4(1:n-1,1);
e5(1,1) = e4(n,1);
del_2 = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1],[-(n-nx) -nx -(nx-1) -1 0 1 nx-1 nx (n-nx)],n,n);
del_2 = del_2/delta_x^2; %relies on delta_x = delta_y

%% Gray-Scott

[t,uv] = ode45(@(t,uv) rhs_GrayScott(uv, del_2, n, ru, rv, f, k), tspan, uv_0,[]);

u = uv(:,1:n);
v = uv(:,n+1:end);

%% Plot Still
figure
subplot(2,1,1)
plot_u = reshape(u(end,:),nx,ny);
contourf(xgrid,ygrid,plot_u);
colorbar
title('u')
subplot(2,1,2)
plot_v = reshape(v(end,:),nx,ny);
contourf(xgrid,ygrid,plot_v);
colorbar
title('v')
return;

%% Plot Movie
figure
% vidw = VideoWriter('gray_scott_movie.avi');
% open(vidw);
movie_frames = [];
for j = 1:size(uv,1)
    frame_u = reshape(u(j,:),nx,ny);
    frame_v = reshape(v(j,:),nx,ny);
%     subplot(2,1,1)
    contourf(xgrid,ygrid,frame_u);
    colorbar
%     subplot(2,1,2)
%     contourf(xgrid,ygrid,frame_v);
%     colorbar
    movie_frames = [movie_frames getframe()];
%     writeVideo(vidw,getframe());
end
% close(vidw);
save gray_scott_movie.mat movie_frames;


