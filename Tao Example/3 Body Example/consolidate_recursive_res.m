all_mr_res = cell(3,1);
all_sep_recon = cell(3,1);
all_params = cell(3,1);
all_r = zeros(3,1);

% Pass 1
load('mwDMD_mr_res_i2.mat'); %use iteration 2 on first pass only
load('mwDMD_sep_recon.mat')
load('mwDMD_params.mat');

sep_recon = struct('nComponents',nComponents,'suppress_growth',suppress_growth,'tspan',tspan,'xr_sep',xr_sep);
params = struct('initialize_artificially',initialize_artificially,'nComponents',nComponents,...
    'nSlide',nSlide,'nSplit',nSplit,'nVars',nVars,'r',r,'stepSize',stepSize,'thresh_pct',...
    thresh_pct,'use_last_freq',use_last_freq,'use_median_freqs',use_median_freqs,...
    'wSteps',wSteps);

all_mr_res{1} = mr_res;
all_sep_recon{1} = sep_recon;
all_params{1} = params;
all_r(1) = r;

% Pass 2
load('mwDMD_mr_res_slow_recursion.mat');
load('mwDMD_sep_recon_slow_recursion.mat')
load('mwDMD_params_slow_recursion.mat');

sep_recon = struct('nComponents',nComponents,'suppress_growth',suppress_growth,'tspan',tspan,'xr_sep',xr_sep);
params = struct('initialize_artificially',initialize_artificially,'nComponents',nComponents,...
    'nSlide',nSlide,'nSplit',nSplit,'nVars',nVars,'r',r,'stepSize',stepSize,'thresh_pct',...
    thresh_pct,'use_last_freq',use_last_freq,'use_median_freqs',use_median_freqs,...
    'wSteps',wSteps);

all_mr_res{2} = mr_res;
all_sep_recon{2} = sep_recon;
all_params{2} = params;
all_r(2) = r;

% Pass 3
load('mwDMD_mr_res_slower_recursion.mat');
load('mwDMD_sep_recon_slower_recursion.mat')
load('mwDMD_params_slower_recursion.mat');

sep_recon = struct('nComponents',nComponents,'suppress_growth',suppress_growth,'tspan',tspan,'xr_sep',xr_sep);
params = struct('initialize_artificially',initialize_artificially,'nComponents',nComponents,...
    'nSlide',nSlide,'nSplit',nSplit,'nVars',nVars,'r',r,'stepSize',stepSize,'thresh_pct',...
    thresh_pct,'use_last_freq',use_last_freq,'use_median_freqs',use_median_freqs,...
    'wSteps',wSteps);

all_mr_res{3} = mr_res;
all_sep_recon{3} = sep_recon;
all_params{3} = params;
all_r(3) = r;


save('recursive_res.mat','all_mr_res','all_sep_recon','all_params','all_r','-v7.3')
