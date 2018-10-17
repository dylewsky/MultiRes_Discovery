close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;

plot_linewidth = 1.3;
slow_linewidth = 2;

load('fig_size_in.mat');
load('example4_mwDMD_sep_recon.mat');
load('example4_raw_data.mat');


figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
nVars = size(x,2);
p_x = cell(nVars,1);
p_xl = cell(nVars,1);
p_xh = cell(nVars,1);
for j = 1:nVars
    if j == 1
        p_x{j} = plot(TimeSpan,x(:,j),'k','LineWidth',plot_linewidth);
        hold on
        p_xl{j} = plot(tspan,xr_L(j,:),'b','LineWidth',slow_linewidth);
        hold on
        p_xh{j} = plot(tspan,xr_H(j,:),'r','LineWidth',plot_linewidth);
        hold on
    else
        p_x{j} = plot(TimeSpan,x(:,j),'k','LineWidth',plot_linewidth,'HandleVisibility','off');
        hold on
        p_xl{j} = plot(tspan,xr_L(j,:),'b','LineWidth',slow_linewidth,'HandleVisibility','off');
        hold on
        p_xh{j} = plot(tspan,xr_H(j,:),'r','LineWidth',plot_linewidth,'HandleVisibility','off');
        hold on
    end
end

set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlabel('Time')
xlim([tspan(1) tspan(end)])
ylim([-6 4])
legend([p_x{1}, p_xl{1}, p_xh{1}],'Input Data','Slow Recon.','Fast Recon.','Location','southeast');
% legend

export_fig '../sep_recon' -pdf -eps -transparent;
