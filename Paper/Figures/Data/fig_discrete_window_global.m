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

load('fig_size_in.mat')
load('example4_raw_data.mat');
load('example4_discrete_window_DMD_mr_res.mat');
load('example4_discrete_window_DMD_params.mat');

wind_edges = wSteps:wSteps:nSteps;
yBounds = [-6 4];
lWidth = 1.4;
inWidth = 2;
reconAlpha = 1;

%% Windowed Recon Figure

figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
rawDataPlots = cell(nVars,1);
for j = 1:nVars
    if j == 1
        rawDataPlots{j} = plot(TimeSpan,x(:,j),'k-','LineWidth',inWidth);
    else
        rawDataPlots{j} = plot(TimeSpan,x(:,j),'k-','LineWidth',inWidth,'HandleVisibility','off');    
    end
    hold on
end
set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlim([TimeSpan(1),TimeSpan(4*wSteps)]);
ylim(yBounds);
xlabel('Time')
hold on
edgePlots = cell(nVars,1);
for nw = 1:length(wind_edges)
    edgePlots{nw} = plot(TimeSpan(wind_edges(nw))*[1 1],yBounds,'k:','HandleVisibility','off');
    hold on
end

reconPlots = cell(nVars, nSlide);
for k = 1:nSlide
    w = mr_res{k}.w;
    b = mr_res{k}.b;
    Omega = mr_res{k}.Omega;
    om_class = mr_res{k}.om_class;
    t = mr_res{k}.t;
    c = mr_res{k}.c;
    t_start = mr_res{k}.t_start;
    tShift = t-t(1); %compute each segment of xr starting at "t = 0"
    xr_window = real(w*diag(b)*exp(Omega*(t-t_start)) + c);
    for j = 1:nVars
        if (j == 1) && (k == 1)
            reconPlots{j,k} = plot(t,xr_window(j,:),'Color',co(2,:),'LineWidth',lWidth);
        else
            reconPlots{j,k} = plot(t,xr_window(j,:),'Color',co(2,:),'LineWidth',lWidth,'HandleVisibility','off');
        end
        hold on
        reconPlots{j,k}.Color(4) = reconAlpha;
    end
    hold on
end
% legend
legend([rawDataPlots{1} reconPlots{1}], 'Input Data', 'Windowed Reconstruction','Location','southeast')
% legend([rawDataPlots{1} reconPlots{1}])
export_fig '../discrete_window_recon' -pdf -eps -transparent;

%% Global Recon Figure
load('example4_mwDMD_sep_recon.mat');
xr_glob = (xr_H + xr_L).';
figure
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
rawDataPlots = cell(nVars,1);
reconPlots = cell(nVars,1);
for j = 1:nVars
    if j == 1
        rawDataPlots{j} = plot(TimeSpan,x(:,j),'k-','LineWidth',inWidth);
        hold on
        reconPlots{j} = plot(tspan,xr_glob(:,j),'Color',co(2,:),'LineWidth',lWidth);
    else
        rawDataPlots{j} = plot(TimeSpan,x(:,j),'k-','LineWidth',inWidth,'HandleVisibility','off');    
        hold on
        reconPlots{j} = plot(tspan,xr_glob(:,j),'Color',co(2,:),'LineWidth',lWidth,'HandleVisibility','off');
    end
    hold on
end
set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlim([TimeSpan(1),TimeSpan(4*wSteps)]);
ylim(yBounds);
xlabel('Time')
hold on

% legend
legend([rawDataPlots{1} reconPlots{1}], 'Input Data', 'Global Reconstruction','Location','southeast')
% legend([rawDataPlots{1} reconPlots{1}])
export_fig '../global_recon' -pdf -eps -transparent;

%% Spectrum Figure

wind_edges = 1:wSteps:nSteps+1;

figure
yBoundsSpec = [0 11];
set(gcf, 'Units', 'Inches', 'Position', fig_size_in)
for nw = 1:length(wind_edges)
    edgePlots{nw} = plot(TimeSpan(wind_edges(nw))*[1 1],yBoundsSpec,'k:','HandleVisibility','off');
    hold on
end
set(gca,'FontSize',fSize)%,'FontWeight','bold')
xlim([TimeSpan(1),TimeSpan(4*wSteps)]);
ylim(yBoundsSpec);
xlabel('Time')
ylabel('Frequency Magnitude |\omega|')

for k = 1:nSlide
    Omega = mr_res{k}.Omega;
    om_class = mr_res{k}.om_class;
    t = mr_res{k}.t;
    wFreqs1 = Omega(om_class == 1);
    wFreqs2 = Omega(om_class == 2);
%     wFreqs1 = wFreqs1 .* conj(wFreqs1);
%     wFreqs2 = wFreqs2 .* conj(wFreqs2);
    wFreqs1 = abs(wFreqs1);
    wFreqs2 = abs(wFreqs2);
    for j = 1:length(wFreqs1)
        plot(TimeSpan(wind_edges([k k+1])), wFreqs1(j) * [1 1], 'b','LineWidth',lWidth);
        hold on
    end
    for j = 1:length(wFreqs2)
        plot(TimeSpan(wind_edges([k k+1])), wFreqs2(j) * [1 1], 'r','LineWidth',lWidth);
        hold on
    end
end
export_fig '../discrete_window_spec' -pdf -eps -transparent;

