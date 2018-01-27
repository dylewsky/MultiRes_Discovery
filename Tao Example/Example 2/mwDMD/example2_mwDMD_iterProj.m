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

slopes_hf = zeros(nSlide,r);
slopes_lf = zeros(nSlide,r);

xr_hf = zeros(r,nSlide);
xr_lf = zeros(r,nSlide);
pcts_hf = zeros(nSlide,1);

for k = 1:nSlide
    w = squeeze(allModes(k,:,:));
    Omega = allFreqs(k,:).';
    hf_modes = w(:,hf_ind);
    lf_modes = w(:,lf_ind);
    hf_modes_orth = GramSchmidt(hf_modes);
    lf_modes_orth = GramSchmidt(lf_modes);
    b = mr_res{k}.b;
    
    mSlopes = zeros(size(w));
    for j = 1:r
        mSlopes(:,j) = w(:,j) * b(j) * Omega(j);
    end
    
    slopes_hf(k,:) = sum(mSlopes(:,hf_ind),2);
    slopes_lf(k,:) = sum(mSlopes(:,lf_ind),2);
    midStep = wMids(k);
    delta_x = [stepSize * delta_t; x(:,midStep+floor(stepSize/2))-x(:,midStep-floor(stepSize/2))].';
    delta_hf = [stepSize * delta_t, slopes_hf(k,:)*stepSize*delta_t];
    delta_lf = [stepSize * delta_t, slopes_lf(k,:)*stepSize*delta_t];
    cosTh_hf = dot(delta_x,delta_hf)/(norm(delta_x)*norm(delta_hf));
    cosTh_lf = dot(delta_x,delta_lf)/(norm(delta_x)*norm(delta_lf));
    pct_hf = cosTh_hf/(cosTh_hf + cosTh_lf);
    pcts_hf(k) = pct_hf;
    if k == 1
        xr_hf(:,k) = x(:,wMids(k)) * pct_hf;
        xr_lf(:,k) = x(:,wMids(k)) * (1-pct_hf);
        continue
    end
    xr_hf(:,k) = xr_hf(:,k-1) + delta_x(2:end).' * pct_hf; %1st element of delta_x is time
    xr_lf(:,k) = xr_lf(:,k-1) + delta_x(2:end).' * (1-pct_hf);
end

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

