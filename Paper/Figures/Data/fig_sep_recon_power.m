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
load('example4_mwDMD_sep_recon.mat');
load('example4_raw_data.mat');

x = sum(x,2);
xr_L = sum(xr_L,1);
xr_H = sum(xr_H,1);

figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
nVars = size(x,2);
N = size(xr_H,2);
Fs = 1/(TimeSpan(2)-TimeSpan(1)); %sampling freq

xdft = fft(x(1:N,:));
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

xdft_H = fft(xr_H.');
xdft_H = xdft_H(1:N/2+1);
psdx_H = (1/(Fs*N)) * abs(xdft_H).^2;
psdx_H(2:end-1) = 2*psdx_H(2:end-1);

xdft_L = fft(xr_L.');
xdft_L = xdft_L(1:N/2+1);
psdx_L = (1/(Fs*N)) * abs(xdft_L).^2;
psdx_L(2:end-1) = 2*psdx_L(2:end-1);

p_ps = plot(freq,10*log10(psdx),'k-','LineWidth',input_linewidth,'DisplayName','Model Data');
hold on
p_psH = plot(freq,10*log10(psdx_H),'r-','LineWidth',plot_linewidth,'DisplayName','Fast Recon.');
hold on
p_psL = plot(freq,10*log10(psdx_L),'b-','LineWidth',plot_linewidth,'DisplayName','Slow Recon.');
% hold on
% p_psLH = plot(freq,10*log10(psdx_L + psdx_H),'g-','LineWidth',plot_linewidth)
grid on
p_psH.Color(4) = recon_alpha;
p_psL.Color(4) = recon_alpha;
% title('Power Spectra of Reconstructions')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
xlim([0 3.5])
legend

set(gca,'FontSize',fSize)%,'FontWeight','bold')
% xlabel('t')
% xlim([tspan(1) tspan(end)])
% ylim([-6 4])
% legend([p_x{1}, p_xl{1}, p_xh{1}],'Input Data','Slow Recon.','Fast Recon.','Location','southeast');
% legend

export_fig '../sep_recon_power' -pdf -eps -transparent;
% print(gcf, '-dpdf', '../sep_recon_power'); 
% print(gcf, '-depsc', '../sep_recon_power'); 
