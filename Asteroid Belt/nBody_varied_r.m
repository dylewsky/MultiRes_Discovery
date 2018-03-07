% clear variables; close all; clc

load('nBody_data.mat');
addpath('optdmd-master');
n = size(z,2)/4;

x = z(:,1:2*n);
% y = z(:,n+1:2*n);
% 
% figure
% plot(t,x(:,1:n))
% hold on
% plot(t,x(:,n+1:2*n))
% hold off
% figure
% plot(t,(x(:,1:n).^2 + x(:,n+1:2*n).^2).^(1/2))
% title('r')
% figure
% plot(t,unwrap(angle(x(:,1:n)+sqrt(-1).*x(:,n+1:2*n))))
% title('\theta')

x = x.';
TimeSpan = t;

% r = 8; %rank to fit w/ optdmd
%  imode = 1, fit full data, slower
%  imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})

% nComponents = 2;

global_SVD = 1; %use global SVD modes for each DMD rather than individual SVD on each window

if global_SVD == 1
    [u,~,v] = svd(x,'econ');
    u = u(:,1:r); %just first r modes
    v = v(:,1:r);
end

initialize_artificially = 0; %this is required for use_last_freq OR use_median_freqs
use_last_freq = 0; 
use_median_freqs = 0; %mutually exclusive w/ use_last_freq

if use_median_freqs == 1
    load('km_centroids.mat');
    freq_meds = repelem(sqrt(km_centroids),r/nComponents);
end

wSteps = 500;
nSplit = floor(length(TimeSpan)/wSteps); %number of windows if they were non-overlapping

% stepSize = wSteps/10;
stepSize = 50;

thresh_pct = 1;

% save('mwDMD_params.mat','r','nComponents','initialize_artificially','use_last_freq','use_median_freqs','wSteps','nSplit','nSteps','nVars','thresh_pct','stepSize','nSlide');

%% execute optDMD

ranks = 4:4:32;

for rn = 3:length(ranks)
    r = ranks(rn);
    nComponents = r/2;
%     figure
    [mr_res, km_centroids] = run_optDMD(x,TimeSpan,r,wSteps,nSplit,stepSize,nComponents,thresh_pct,use_last_freq,use_median_freqs);
    plot_optDMD(mr_res,km_centroids,x,TimeSpan,r,wSteps,nSplit,stepSize,nComponents,thresh_pct);
    outFile = ['varied_r_' num2str(rn) '.fig'];
    save(['varied_r_res_' num2str(rn) '.mat'],'mr_res','km_centroids');
    saveas(gcf,outFile);
end

