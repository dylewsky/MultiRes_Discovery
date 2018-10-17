close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;


centroid_darkness = 0.9;
% centroid_linewidth = 5;
% centroid_alpha = 0.6;
centroid_linewidth = 2;
centroid_alpha = 1;

load('fig_size_in.mat')
load('Three_Body_om_spec.mat');
load('Three_Body_Data_Cartesian.mat');

all_om = reshape(real(omega_series),numel(omega_series),1);
all_om = all_om(isnan(all_om) == 0);
[~, km_centroids] = kmeans(all_om,3,'Distance','cityblock','Replicates',5);
km_centroids = sort(km_centroids);

vert_lims = [0 2.3*10^4];
% centroid_lims = [0 3400];
centroid_lims = vert_lims;

figure
omHist = histogram(all_om,48,'FaceColor',0.1*[1 1 1]);
hold on
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
cent1 = plot([km_centroids(1) km_centroids(1)],centroid_lims,'Color',[0 0 1]*centroid_darkness,'LineWidth',centroid_linewidth);
hold on
cent2 = plot([km_centroids(2) km_centroids(2)],centroid_lims,'Color',[1 0 0]*centroid_darkness,'LineWidth',centroid_linewidth);
hold on
cent3 = plot([km_centroids(3) km_centroids(3)],centroid_lims,'Color',[0 0.7 0]*centroid_darkness,'LineWidth',centroid_linewidth);
cent1.Color(4) = centroid_alpha;
cent2.Color(4) = centroid_alpha;
cent3.Color(4) = centroid_alpha;
hold off

uistack(omHist,'top'); %draw histogram on top of centroids

set(gca,'FontSize',fSize)%,'FontWeight','bold')
% xlim([-5 125])
ylim(vert_lims)
xlabel('|\omega_i^k|^2','Interpreter','tex','FontWeight','bold')
ylabel('Count')

labelSize = 12;
% c1_label = text(20,7000,{'Cluster #1', 'Centroid'},'FontSize',labelSize,'HorizontalAlignment','center',...
%     'Color',[0 0 1]*centroid_darkness);
c1_label = annotation('textarrow',0.18+[0.12 0],0.6+[0.07 0],'String',{'Cluster #1', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[0 0 1]*centroid_darkness,...
    'Color',[0 0 1]*centroid_darkness);

c2_label = annotation('textarrow',0.29+[0.07 0],0.84+[0 0.03],'String',{'Cluster #2', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[1 0 0]*centroid_darkness,...
    'Color',[1 0 0]*centroid_darkness);

c3_label = annotation('textarrow',0.78+[0 0.07],0.84+[0 0.03],'String',{'Cluster #3', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[0 0.7 0]*centroid_darkness,...
    'Color',[0 0.7 0]*centroid_darkness);

% figure
% semilogy(time_series,omega_series(1:2,:),'b.')
% hold on
% semilogy(time_series,omega_series(3:4,:),'r.')
% hold off

set(gcf,'color','w');
export_fig '../om_hist_3body' -pdf -eps -transparent;

% print(gcf, '-dpdf', '../om_hist_3body'); 
% print(gcf, '-depsc', '../om_hist_3body'); 