clear variables; close all; clc

load('../raw_data_2_hiRes.mat');
load('mwDMD_mr_res.mat');
load('mwDMD_allModes.mat');
load('mwDMD_params.mat');

%%
wMids = stepSize*(0:nSlide-1) + ceil(wSteps/2);
delta_t = TimeSpan(2)-TimeSpan(1);
tMids = wMids * delta_t;

lf_ind = [1 2];
hf_ind = [3 4];

x_hf = zeros(1,nSlide);
x_lf = zeros(1,nSlide);

for k = 1:nSlide
    xW = mr_res{k}.x;
    xMean = mean(xW,2);
    w = squeeze(allModes(k,:,:));
    Omega = allFreqs(k,:).';
    hf_modes = w(:,hf_ind);
    lf_modes = w(:,lf_ind);
    b = mr_res{k}.b;
    if (norm(hf_modes(:,1)-conj(hf_modes(:,2))) >= 10^(-6)) || (norm(lf_modes(:,1)-conj(lf_modes(:,2))) >= 10^(-6))
        disp('not conjugate pair')
        continue
    end
    hf_mode = real(hf_modes(:,1));
    lf_mode = real(lf_modes(:,1));
    xdh = dot(hf_mode,xMean);
    xdl = dot(lf_mode,xMean);
    hdh = dot(hf_mode,hf_mode);
    ldl = dot(lf_mode,lf_mode);
    hdl = dot(lf_mode,hf_mode);
    x_hf(k) = (hdl*xdl - xdh*ldl)/(hdl^2 - hdh*ldl); %non-orthogonal projection
    x_lf(k) = (hdl*xdh - xdl*hdh)/(hdl^2 - hdh*ldl);
end

%%
figure
subplot(2,1,1)
plot(TimeSpan(1:nSteps),x(:,1:nSteps));
subplot(2,1,2)
plot(tMids,x_hf,'r-',tMids,x_lf,'b-');

%%
figure('units','pixels','Position',[100 100 1366 768])

tk = delta_t * (-floor(stepSize/2):floor(stepSize/2));
for j = 1:r
    subplot(2,2,j)
    plot(TimeSpan(1:nSteps),x(j,1:nSteps),'k-')
    hold on
    for k = 1:nSlide
        midStep = wMids(k);
        plot(tk+TimeSpan(midStep),real(slopes_hf(k,j))*tk + x(j,midStep),'r-')
        hold on
        plot(tk+TimeSpan(midStep),real(slopes_lf(k,j))*tk + x(j,midStep),'b-')
        hold on
    end
end

%%
pcts_col = zeros(1,nSteps);
n_col = zeros(1,nSteps);
for k = 1:nSlide
    pcts_col(wMids(k)-floor(wSteps/2)+1:wMids(k)+floor(wSteps/2)) = pcts_col(wMids(k)-floor(wSteps/2)+1:wMids(k)+floor(wSteps/2)) + pcts_hf(k);
    n_col(wMids(k)-floor(wSteps/2)+1:wMids(k)+floor(wSteps/2)) = n_col(wMids(k)-floor(wSteps/2)+1:wMids(k)+floor(wSteps/2)) + 1;
end
pcts_col = real(pcts_col./n_col);
% pcts_col = [zeros((nSteps - nSlide*stepSize)/2,1); repelem(real(pcts_hf),stepSize); zeros((nSteps - nSlide*stepSize)/2,1)].';
figure('units','pixels','Position',[100 100 1366 768])

subplot(3,1,1)
for j = 1:r
    surface([TimeSpan(1:nSteps);TimeSpan(1:nSteps)],[x(j,1:nSteps);x(j,1:nSteps)],zeros(2,nSteps),[pcts_col; pcts_col],...
            'facecol','no','edgecol','interp','linew',2);
    hold on
end
cm = colormap('winter');
caxis([0 1]);
colorbar
subplot(3,1,2)
plot(tMids,real(xr_hf),'Color',cm(end,:))
title('HF Recon.')
subplot(3,1,3)
plot(tMids,real(xr_lf),'Color',cm(1,:))
title('LF Recon.')

