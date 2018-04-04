close all; clear variables; clc

addpath('../../../diffusion_maps-master');

load('modeSeries_i2.mat');

%% Diffusion Map on Full Mode Data

X = [real(modeStack(:,[1:4,9:12])), imag(modeStack(:,[1:4,9:12]))];

N = size(X,1);
K = zeros(N);

alphas = 3; %determines kernel scale

for n = 1:length(alphas)
    alpha = alphas(n);
    for j = 1:N
        for k = 1:j
            pxy = exp(-norm(X(j,:)-X(k,:))^2/alpha);
            K(j,k) = pxy;
            K(k,j) = pxy;
        end
    end

    [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(K);
    
end

save('diffusion_map_full.mat','diffusion_map');

%% Diffusion Map on Separated Mode Data

X_hf = [real(modeStack(:,1:4)), imag(modeStack(:,1:4))];
X_lf = [real(modeStack(:,9:12)), imag(modeStack(:,9:12))];

N = size(X_hf,1);
K_hf = zeros(N);
K_lf = zeros(N);

alphas = 3; %determines kernel scale

for n = 1:length(alphas)
    alpha = alphas(n);
    for j = 1:N
        for k = 1:j
            pxy_hf = exp(-norm(X_hf(j,:)-X_hf(k,:))^2/alpha);
            pxy_lf = exp(-norm(X_lf(j,:)-X_lf(k,:))^2/alpha);
            K_hf(j,k) = pxy_hf;
            K_hf(k,j) = pxy_hf;
            K_lf(j,k) = pxy_lf;
            K_lf(k,j) = pxy_lf;
        end
    end

    [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(K_hf);
    title('HF')
    save('diffusion_map_hf.mat','diffusion_map');
    
    [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(K_lf);
    title('LF')
    save('diffusion_map_lf.mat','diffusion_map');
end
