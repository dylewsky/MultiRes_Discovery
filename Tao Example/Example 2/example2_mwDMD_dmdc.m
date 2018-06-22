clear variables; close all; clc

load mwDMD_sep_recon.mat;
load mwDMD_params.mat;
load '../raw_data_2_hiRes.mat';

r = 4;

xr_H = real(xr_H(:,1:end-stepSize));
xr_L = real(xr_L(:,1:end-stepSize));

dmdc_res = cell(2,1);

for cVar = 1:2
    if cVar == 1
        StateData = xr_H;
        InputData = xr_L;
    else
        StateData = xr_L;
        InputData = xr_H;
    end
    X = StateData(:,1:end-1);
    Xp = StateData(:,2:end);
    Ups = InputData(:,1:end-1);
    
    Omega = [X; Ups];
    
    [Util, Sigtil, Vtil] = svd(Omega,'econ');
    [Uhat, Sighat, Vbar] = svd(Xp,'econ');
    
    n = size(X,1);
    q = size(Ups,1);
    
    U_1 = Util(1:n,:);
    U_2 = Util(n+1:n+q,:);
    
    approxA = Uhat' * Xp * Vtil * inv(Sigtil) * U_1' * Uhat;
    approxB = Uhat' * Xp * Vtil * inv(Sigtil) * U_2';
    
    [W,D] = eig(approxA);
    Phi = Xp * Vtil * inv(Sigtil) * U_1' * Uhat * W;
    
    dmdc_res{cVar}.Phi = Phi;
    dmdc_res{cVar}.approxA = approxA;
    dmdc_res{cVar}.approxB = approxB;
    figure('Position',[200 200 1000 400])
    subplot(2,3,1:2)
    plot(TimeSpan(1:size(xr_H,2)),StateData)
    title('State Data')
    annotation('textbox',[.51 .61 .12 .05],'String','x_{k+1} = Ax_k + Bu_k','LineStyle','none')
    subplot(2,3,4:5)
    plot(TimeSpan(1:size(xr_H,2)),InputData)
    title('Control Input')
    subplot(2,3,3)
    imagesc(abs(approxA - eye(r)))
    colorbar
    title('|A - I|')
    subplot(2,3,6)
    imagesc(abs(approxB))
    colorbar
    title('|B|')
    
    % reconstruct from results
    xr = zeros(size(StateData));
    xr(:,1) = StateData(:,1);
    for tn = 2:size(StateData,2)
        xr(:,tn) = approxA * xr(:,tn-1) + approxB * InputData(:,tn-1);
    end
    figure
    plot(TimeSpan(1:size(StateData,2)),StateData,'k');
    hold on
    plot(TimeSpan(1:size(StateData,2)),xr,'r');
    hold off
    title('State Data vs. Reconstruction')
    ylim(1.5*[min(min(StateData)) max(max(StateData))]);
    
    dmdc_res{cVar}.xr = xr;
end