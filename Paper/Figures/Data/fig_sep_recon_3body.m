close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;
comp_colors = cell(nComponents,1);
comp_colors{1} = [0 0 1];
comp_colors{2} = [1 0 0];
comp_colors{3} = [0 0.7 0];

plot_linewidth = 2*[1 1 1 1.5]; %comp. 1, 2 , 3, raw data

load('fig_size_in.mat');
load('Three_Body_Data_Cartesian.mat');
TimeSpan = tspan;
load('Three_Body_mwDMD_sep_recon.mat');
load('Three_Body_mwDMD_params.mat');


figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
p_x = cell(nVars,1);
p_xr = cell(nVars,nComponents);
% p_x = cell(nVars,1);
% p_xl = cell(nVars,1);
% p_xh = cell(nVars,1);
for j = 1:nVars
    if j == 1
        p_x{j} = plot(TimeSpan,pos(j,:),'k','LineWidth',plot_linewidth(4));
        hold on
        for k = 1:nComponents
                p_xr{j,k} = plot(tspan,xr_sep{k}(j,:),'Color',comp_colors{k},'LineWidth',plot_linewidth(k));
                hold on
        end
    else
        p_x{j} = plot(TimeSpan,pos(j,:),'k','LineWidth',plot_linewidth(4),'HandleVisibility','off');
        hold on
        for k = 1:nComponents
                p_xr{j,k} = plot(tspan,xr_sep{k}(j,:),'Color',comp_colors{k},'LineWidth',plot_linewidth(k),'HandleVisibility','off');
                hold on
        end
    end
end

set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlabel('t (Years)')
ylabel('Cartesian Planet Positions (AU)')
% xlim([tspan(1) tspan(end)])
xlim([10 10.15]*10^3);
ylim([-20 11])
legend([p_x{1}, p_xr{1,1}, p_xr{1,2}, p_xr{1,3}],'Input Data','Component 1: Slow','Component 2: Saturn','Component 3: Jupiter','Location','southeast');
% legend

export_fig '../sep_recon_3body' -pdf -eps -transparent;
