close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;
comp_colors = cell(nComponents,1);
comp_colors{1} = [0 0 1];
comp_colors{2} = [1 0 0];
comp_colors{3} = [0 0.7 0];
co = [0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0    0.4470    0.7410;
    0.6500    0.3250    0.0980;
    0.4660    0.6740    0.1880;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

plot_linewidth = 2*[1 1 1 1.5]; %comp. 1, 2 , 3, raw data

load('fig_size_in.mat');
load('Three_Body_Data_Cartesian.mat');
load('Three_Body_Data_Orbital.mat');
TimeSpan = tspan;
load('Three_Body_mwDMD_sep_recon.mat');
load('Three_Body_mwDMD_params.mat');


figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
slow_1 = [OOmega(1,:); aa(1,:); ee(1,:); ii(1,:); oomega(1,:)];
slow_2 = [OOmega(2,:); aa(2,:); ee(2,:); ii(2,:); oomega(2,:)];
for j = 1:5 %normalize
    slow_1(j,:) = slow_1(j,:)/norm(slow_1(j,:));
    slow_1(j,:) = slow_1(j,:)-mean(slow_1(j,:));
    slow_2(j,:) = slow_2(j,:)/norm(slow_2(j,:));
    slow_2(j,:) = slow_2(j,:)-mean(slow_2(j,:));
end
plot(tcoarse,slow_1,'k')
hold on
plot(tcoarse,slow_2,'b')
hold off

set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlabel('t (Years)')
ylabel('Cartesian Planet Positions (AU)')
% xlim([tspan(1) tspan(end)])
xlim([tcoarse(1) tcoarse(end)]);
% ylim([-20 11])
% legend([p_x{1}, p_xr{1,1}, p_xr{1,2}, p_xr{1,3}],'Input Data','Component 1: Slow','Component 2: Saturn','Component 3: Jupiter','Location','southeast');
% legend

slow_all = [slow_1;slow_2];
[Uk,Sk,Vk] = svd(slow_all,'econ');

slow_dmd = xr_sep{1};
for j = 1:nVars %normalize
    slow_dmd(j,:) = slow_dmd(j,:)/norm(slow_dmd(j,:));
    slow_dmd(j,:) = slow_dmd(j,:)-mean(slow_dmd(j,:));
end
[Ud,Sd,Vd] = svd(slow_dmd,'econ');

figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
plot(cumsum(diag(Sk))/sum(diag(Sk)),'k.-','LineWidth',2,'MarkerSize',20)
% plot(diag(Sk),'ko-','LineWidth',2)
hold on
plot(cumsum(diag(Sd))/sum(diag(Sd)),'.-','Color',co(2,:),'LineWidth',2,'MarkerSize',20)
% plot(diag(Sd),'ro-','LineWidth',2)
legend('Slow Keplerian Orbital Elements (e, a, i, \Omega, \omega)', ...
    'Slow-Scale DMD Reconstruction','Location','southeast');
xlim([1 length(diag(Sk))])
set(gca,'XTick',1:length(diag(Sk)))
ylim([0 1]);
xlabel('Mode #')
ylabel('Singular Value (Normalized Cumulative Sum)')
set(gca,'FontSize',fSize)%,'FontWeight','bold')
export_fig '../3body_orbital_param_SV' -pdf -eps -transparent;


SVk = Sk*Vk.';
SVd = Sd*Vd.';
SVk = SVk(1:6,:);
tcoarse = tcoarse(1:99590);
SVk = SVk(:,1:99590);
SVd_interp = interp1(tspan,SVd.',tcoarse).';

% SVd_interp = SVd_interp(:,1:99590);

figure
for j = 1:nVars
    if(norm(SVk(j,:) - SVd_interp(j,:)) > norm(SVk(j,:) + SVd_interp(j,:)))
        SVk(j,:) = -SVk(j,:); %sync signs of SVD results
        disp(['Flip mode ' num2str(j)])
    end
    subplot(nVars,1,j)
    plot(tcoarse,SVk(j,:),'k')
    hold on
    plot(tcoarse,SVd_interp(j,:),'b')
end


figure
plot(tcoarse,SVk.','k')
hold on
plot(tcoarse,SVd_interp.','b','LineWidth',2)
hold off



% export_fig '../3body_orbital_param' -pdf -eps -transparent;
