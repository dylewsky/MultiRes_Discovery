close all;

addpath('altmany-export_fig-9ac0917');
% addpath(genpath('CUSTOMPATH/cvx'))
% rmpath('CUSTOMPATH/cvx/lib/narginchk_:')

% Plot Style
fSize = 12;

plot_linewidth = 2;
input_linewidth = 4;

recon_alpha = 1;

load('fig_size_in.mat');
load('Three_Body_mwDMD_sep_recon.mat');
load('Three_Body_Data_Cartesian.mat');
load('Three_Body_mwDMD_params.mat');
TimeSpan = tspan;

x = sum(pos,1).';
xr = cell(size(xr_sep));
for j = 1:nComponents
    xr{j} = sum(xr_sep{j},1).';
end

figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
N = size(xr{1},1);
Fs = 1/(TimeSpan(2)-TimeSpan(1)); %sampling freq

xdft = fft(x(1:N,:));
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

p_ps = plot(freq,10*log10(psdx),'k-','LineWidth',input_linewidth,'DisplayName','Model Data');
hold on

xdft_R = cell(size(xr));
psdx_R = cell(size(xr));
p_psR = cell(size(xr));

displayNames = {'Slow','Saturn','Jupiter'};
colors = {[0 0 1],[1 0 0],[0 0.7 0]};

for j = 1:nComponents
    xdft_R{j} = fft(xr{j});
    xdft_R{j} = xdft_R{j}(1:N/2+1);
    psdx_R{j} = (1/(Fs*N)) * abs(xdft_R{j}).^2;
    psdx_R{j}(2:end-1) = 2*psdx_R{j}(2:end-1);
    p_psR{j} = plot(freq,10*log10(psdx_R{j}),'Color',colors{j},'LineWidth',plot_linewidth,'DisplayName',displayNames{j});
    hold on
    p_psR{j}.Color(4) = recon_alpha;
end
grid on
% title('Power Spectra of Reconstructions')
xlabel('Frequency (Years^{-1})')
ylabel('Power/Frequency (dB Years)')
xlim([0 0.15])
legend

set(gca,'FontSize',fSize)%,'FontWeight','bold')
% xlabel('t')
% xlim([tspan(1) tspan(end)])
% ylim([-6 4])
% legend([p_x{1}, p_xl{1}, p_xh{1}],'Input Data','Slow Recon.','Fast Recon.','Location','southeast');
% legend

export_fig '../sep_recon_power_3body' -pdf -eps -transparent;
% print(gcf, '-dpdf', '../sep_recon_power'); 
% print(gcf, '-depsc', '../sep_recon_power'); 
