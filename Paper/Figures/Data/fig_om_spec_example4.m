close all;

addpath('altmany-export_fig-9ac0917');

% Plot Style
fSize = 12;


mSize = 5;
centroid_darkness = 0.4;
centroid_linewidth = 3;
centroid_alpha = 0.8;

load('fig_size_in.mat')
load('example4_om_spec.mat');
load('example4_raw_data.mat');

figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)

all_om = reshape(real(omega_series),numel(omega_series),1);
all_om = all_om(isnan(all_om) == 0);
[~, km_centroids] = kmeans(all_om,2,'Distance','cityblock','Replicates',5);
km_centroids = sort(km_centroids);
cent1 = plot([TimeSpan(1)-2 TimeSpan(end)+1], [km_centroids(1) km_centroids(1)],'Color',[0 0 1]*centroid_darkness,'LineWidth',centroid_linewidth);
hold on
cent2 = plot([TimeSpan(1)-2 TimeSpan(end)+1], [km_centroids(2) km_centroids(2)],'Color',[1 0 0]*centroid_darkness,'LineWidth',centroid_linewidth);
hold on
cent1.Color(4) = centroid_alpha;
cent2.Color(4) = centroid_alpha;

plot(time_series,omega_series(1:2,:),'b.','MarkerSize',mSize)
hold on
plot(time_series,omega_series(3:4,:),'r.','MarkerSize',mSize)
hold off
xlim([TimeSpan(1)-2 TimeSpan(end)+20]);

set(gca,'FontSize',fSize)%,'FontWeight','bold')
ylim([-10 140])
ylabel('|\omega_i^k|^2','Interpreter','tex','FontWeight','bold')
xlabel('Time (Center of DMD Window)')

labelSize = 12;
label_xpos = 0.70+[0.03 0];
c1_label = annotation('textarrow',label_xpos,0.21+[0.05 0],'String',{'Cluster #1', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[0 0 1]*centroid_darkness,...
    'Color',[0 0 1]*centroid_darkness);

c2_label = annotation('textarrow',label_xpos,0.66+[0 0.05],'String',{'Cluster #2', 'Centroid'},...
    'FontSize',labelSize,'HorizontalAlignment','center','TextColor',[1 0 0]*centroid_darkness,...
    'Color',[1 0 0]*centroid_darkness);

% figure
% semilogy(time_series,omega_series(1:2,:),'b.')
% hold on
% semilogy(time_series,omega_series(3:4,:),'r.')
% hold off

export_fig '../om_spec' -pdf -eps -transparent;
