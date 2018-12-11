load('mwDMD_sep_recon_recursion_01.mat')

ds_factor = 100;

pos = xr_sep{1}(:,1:ds_factor:end);
tspan = tspan(1:ds_factor:end);

save('mwDMD_sep_recon_recursion_01_downsample.mat','pos','tspan');