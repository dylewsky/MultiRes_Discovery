close all

f=.055;
k=.062;

epsilon = 0.1;
gamma = 1;

sigma = 1;
lambda = 1;
kappa = 1;

% Diffusion rates
da = 1;
db = .5;
du = 1;
dv = 1;
% Size of grid
width = 128;
% 5,000 simulation seconds with 4 steps per simulated second
dt = .25;
stoptime = 1000;
  
[t, A, B, U, V] = initial_conditions(width);

uMeans = [];
uMeans = [uMeans mean(mean(U))];

lB = my_dx(-my_dx(B)) - my_dy(my_dy(B));

hi = image(U);
hi.CDataMapping = 'scaled';
colorbar
ht = text(3,width-3,'Time = 0');
ht.Color = [.95 .2 .8];

tic
nframes = 1;
while t<stoptime
    anew = A + (da*my_laplacian(A) - A.*B.^2 + f*(1-A))*dt;% + epsilon*(U-V)*dt;
    bnew = B + (db*my_laplacian(B) + A.*B.^2 - (k+f)*B)*dt;% + epsilon*(-U+V)*dt;
%     unew = U + (du*my_laplacian(U) + lambda * U - U.^3/3 - kappa - sigma * V)*dt; % + gamma*(A-B).^2*dt;
%     vnew = V + (dv*my_laplacian(V) + U - V)*dt;
    unew = U + (- sqrt(-1)*U + V)*dt + (gamma*my_d2x(A).^2)*dt;
    vnew = V + (- sqrt(-1)*V - U)*dt + (gamma*my_d2x(B).^2)*dt;
    
    
    A = anew;
    B = bnew;
    lB = my_dx(-my_dx(B)) - my_dy(my_dy(B));

    U = unew;
    V = vnew;
    
%     uMeans = [uMeans mean(mean(U))];
    
    t = t+dt;
%     hi.CData = B;
%     ht.String = ['Time = ' num2str(t)];
%     drawnow limitrate
    hi.CData = lB;
    ht.String = ['Time = ' num2str(t)];
    drawnow limitrate
    nframes = nframes+1;
end

% figure
% plot(uMeans)

return

axes('Position',[0 0 1 1])
axis off

% Add a scaled-color image
hi = image(B);
hi.CDataMapping = 'scaled';


delta = toc;
disp([num2str(nframes) ' frames in ' num2str(delta) ' seconds']);