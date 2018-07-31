close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;


centroid_darkness = 0.8;
centroid_linewidth = 5;
centroid_alpha = 0.6;

load('example4_om_spec.mat');
load('example4_raw_data.mat');

all_om = reshape(real(omega_series),numel(omega_series),1);
all_om = all_om(isnan(all_om) == 0);
[~, km_centroids] = kmeans(all_om,2,'Distance','cityblock','Replicates',5);
km_centroids = sort(km_centroids);

vert_lims = [0 8000];
% centroid_lims = [0 3400];
centroid_lims = vert_lims;

figure
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 5, 4])
histogram(all_om,128,'FaceColor',0.1*[1 1 1])
hold on
cent1 = plot([km_centroids(1) km_centroids(1)],centroid_lims,'Color',[0 0 1]*centroid_darkness,'LineWidth',centroid_linewidth);
hold on
cent2 = plot([km_centroids(2) km_centroids(2)],centroid_lims,'Color',[1 0 0]*centroid_darkness,'LineWidth',centroid_linewidth);
hold off
cent1.Color(4) = centroid_alpha;
cent2.Color(4) = centroid_alpha;
set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlim([-5 125])
ylim(vert_lims)
xlabel('|\omega_i|^2','Interpreter','tex','FontWeight','bold')
ylabel('Count')

labelSize = 12;
% c1_label = text(20,7000,{'Cluster #1', 'Centroid'},'FontSize',labelSize,'HorizontalAlignment','center',...
%     'Color',[0 0 1]*centroid_darkness);
c1_label = annotation('textarrow',0.22+[0.05 0],[0.77 0.70],'String',{'Cluster #1', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[0 0 1]*centroid_darkness,...
    'Color',[0 0 1]*centroid_darkness);

c2_label = annotation('textarrow',0.7+[0 0.05],[0.77 0.70],'String',{'Cluster #2', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[1 0 0]*centroid_darkness,...
    'Color',[1 0 0]*centroid_darkness);

% figure
% semilogy(time_series,omega_series(1:2,:),'b.')
% hold on
% semilogy(time_series,omega_series(3:4,:),'r.')
% hold off

export_fig '../om_hist' -pdf -eps;
