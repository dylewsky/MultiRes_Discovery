close all;

addpath('altmany-export_fig-9ac0917');
addpath('PeterRochford-MarkerTransparency-6a28bba');

% Plot Style
fSize = 12;

co = [0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0    0.4470    0.7410;
    0.6500    0.3250    0.0980;
    0.4660    0.6740    0.1880;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co)

% Figure
load('example4_raw_data.mat');
load('example4_mwDMD_params.mat');
figure
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 27/8, 3], 'PaperUnits', 'Inches', 'PaperSize', [8.5 11])
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 5, 4])
plot([TimeSpan(1) TimeSpan(end)],[0 0],'k-')
hold on
plot(TimeSpan,x,'LineWidth',1.5);
xlim([TimeSpan(1) TimeSpan(end)]);
xlabel('Time (s)')
% set(gca,'YTickLabel','')
% title('Simple Multiscale Model')
set(gca,'FontSize',fSize)%,'FontWeight','bold')
ylim([-7 4]);

ax = gca;
% axis tight

% annotation('line',0.1 + [0 wSteps-1]*(TimeSpan(2)-TimeSpan(1))/(TimeSpan(end)-TimeSpan(1)))
hold on
xStart = 8; %horizontal position to start the bracket
windMultiplier = 2; %factor to increase wSteps by for visualization purposes
bracketX = xStart + [0 windMultiplier*wSteps-1]*(TimeSpan(2)-TimeSpan(1));
bracketY = [-4.5 4];
bracketHeight = 0.2;
ghostTrailDist = diff(bracketX)/3;
ghostTrailAlphas = [0.3 0.15 0.05];
bracketBottom = plot(bracketX, [bracketY(1) bracketY(1)],'k.-','LineWidth',2);
hold on
plot([bracketX(1) bracketX(1)],[bracketY(1) bracketY(1)+bracketHeight],'Color',bracketBottom.Color,'LineStyle',bracketBottom.LineStyle,'LineWidth',bracketBottom.LineWidth)
hold on
plot([bracketX(2) bracketX(2)],[bracketY(1) bracketY(1)+bracketHeight],'Color',bracketBottom.Color,'LineStyle',bracketBottom.LineStyle,'LineWidth',bracketBottom.LineWidth)
% annotation('rectangle', [bracketX(1)/(TimeSpan(end)-TimeSpan(1)), 0.5/9, bracketX(2)/(TimeSpan(end)-TimeSpan(1)), 1],'FaceColor','blue','FaceAlpha',.2)
hold on
wRect = rectangle('Position',[bracketX(1), bracketY(1), bracketX(2)-bracketX(1),bracketY(2)-bracketY(1)],'FaceColor','k','LineStyle','none');
wRect.FaceColor(4) = 0.1;

ghostBottoms = cell(length(ghostTrailAlphas)+1,1);
ghostLefts = ghostBottoms;
ghostRights = ghostBottoms;
ghostBottoms{1} = bracketBottom;
for gh = 1:length(ghostTrailAlphas)
    ghostBottoms{gh+1} = plot(bracketX-gh*ghostTrailDist, [bracketY(1) bracketY(1)],'k.-','LineWidth',ghostBottoms{gh}.LineWidth,'MarkerSize',1);
    ghostBottoms{gh+1}.Color(4) = ghostTrailAlphas(gh);
    ghostBottoms{gh+1}.LineWidth = 1;
    ghostBottoms{gh+1}.MarkerFaceColor = (1 - 1.5*ghostTrailAlphas(gh)) * [1 1 1];
    ghostBottoms{gh+1}.MarkerEdgeColor = (1 - 1.5*ghostTrailAlphas(gh)) * [1 1 1];
    hold on
    ghostLefts{gh+1} = plot([bracketX(1)-gh*ghostTrailDist bracketX(1)-gh*ghostTrailDist],[bracketY(1) bracketY(1)+bracketHeight],'Color',ghostBottoms{gh+1}.Color,'LineStyle',ghostBottoms{gh+1}.LineStyle,'LineWidth',ghostBottoms{gh+1}.LineWidth);
    hold on
    ghostRights{gh+1} = plot([bracketX(2)-gh*ghostTrailDist bracketX(2)-gh*ghostTrailDist],[bracketY(1) bracketY(1)+bracketHeight],'Color',ghostBottoms{gh+1}.Color,'LineStyle',ghostBottoms{gh+1}.LineStyle,'LineWidth',ghostBottoms{gh+1}.LineWidth);
    ghostLefts{gh+1}.Color(4) = ghostTrailAlphas(gh);
    ghostRights{gh+1}.Color(4) = ghostTrailAlphas(gh);
    setMarkerColor(ghostBottoms{gh+1},'k',ghostTrailAlphas(gh));
    hold on
end
for gh = 1:length(ghostTrailAlphas) %move to top in reverse order
    uistack(ghostBottoms{length(ghostTrailAlphas)-gh+1},'top');
end

annotation('arrow',[0.37 0.5],[0.32 0.32],'LineWidth',2,'HeadStyle','vback2')
annotation('textbox',[0.15 0.14 0.2 0.15],...
    'String','$\approx \sum_i^r \vec{\phi}_i e^{\omega_i t} b_i$','Interpreter','latex',...
    'FontSize',16,'LineStyle','none','VerticalAlignment','bottom','FontWeight','bold')
% curlyBrace = text(mean(bracketX),-5,'{','FontSize',36,'HorizontalAlignment','center','VerticalAlignment','middle');
%     'String','asdf','Interpreter','tex','FontSize',16,'LineStyle','none','VerticalAlignment','bottom','FontWeight','bold');
% set(curlyBrace,'Rotation',90)
curlyBrace = drawbrace([2 -5.6],[18 -5.6]);
set(curlyBrace,'Color','k','LineWidth',1)
export_fig '../moving_window' -pdf -eps;
