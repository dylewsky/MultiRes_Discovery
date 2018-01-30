clear variables;
load('res_list_2.mat');
load('mr_res_2.mat');
load('raw_data_2.mat');


sepTrial = 17; %index of trial with window size that successfully separated scales
lfScale = 1;
hfScale = 2; % 1 = low-frequency, 2 = high-frequency
% sepModes = 2; % # of modes 

nVars = size(x,1);
nSteps = res_list(sepTrial,3) * res_list(sepTrial,4); % # windows * steps per window

j = res_list(sepTrial,2);
pn = res_list(sepTrial,1);
nSplit = 2^(j-1);

hf_res = cell(nSplit, 1);
lf_res = cell(nSplit, 1);
for k = 1:nSplit
    om_class = mr_res{pn,j,k}.om_class;
    hf_res{k}.t = mr_res{pn,j,k}.t;
    hf_res{k}.x = mr_res{pn,j,k}.x;
    hf_res{k}.w = mr_res{pn,j,k}.w(:,om_class == hfScale);
    hf_res{k}.Omega = mr_res{pn,j,k}.Omega(om_class == hfScale);
    hf_res{k}.b = mr_res{pn,j,k}.b(om_class == hfScale);
    hf_res{k}.om_post = mr_res{pn,j,k}.om_post(om_class == hfScale,:);
    
    lf_res{k}.t = mr_res{pn,j,k}.t;
    lf_res{k}.x = mr_res{pn,j,k}.x;
    lf_res{k}.w = mr_res{pn,j,k}.w(:,om_class == lfScale);
    lf_res{k}.Omega = mr_res{pn,j,k}.Omega(om_class == lfScale);
    lf_res{k}.b = mr_res{pn,j,k}.b(om_class == lfScale);
    lf_res{k}.om_post = mr_res{pn,j,k}.om_post(om_class == lfScale,:);
end

%% Plot Separated Reconstructions
close all;
figure

t_full = [];
x_full = [];
xr_full = [];
xr_hf_full = [];
xr_lf_full = [];


for k = 1:nSplit
    t = mr_res{pn,j,k}.t;
    x = mr_res{pn,j,k}.x;
    t_full = [t_full t];
    x_full = [x_full x];
    
    w = mr_res{pn,j,k}.w;
    b = mr_res{pn,j,k}.b;
    Omega = mr_res{pn,j,k}.Omega;
    
    bt_window = diag(b)*exp(Omega*t); %time series recon. in mode space
    xr_window = w*bt_window; %time series recon. in original variables
    xr_full = [xr_full xr_window];
    
    w = hf_res{k}.w;
    b = hf_res{k}.b;
    Omega = hf_res{k}.Omega;
    xr_hf = w*diag(b)*exp(Omega*t);
    xr_hf_full = [xr_hf_full xr_hf];
    
    w = lf_res{k}.w;
    b = lf_res{k}.b;
    Omega = lf_res{k}.Omega;
    xr_lf = w*diag(b)*exp(Omega*t);
    xr_lf_full = [xr_lf_full xr_lf];
end

subplot(2,1,1)
pRaw = cell(nVars,1);
pRecon = cell(nVars,1);
for ip = 1:nVars
    pRaw{ip} = plot(t_full,x_full(ip,:),'k','LineWidth',1);
    hold on
    pRecon{ip} = plot(t_full,real(xr_full(ip,:)),'g','LineWidth',1);
end
subplot(2,1,2)
pHF = cell(nVars,1);
pLF = cell(nVars,1);
for ip = 1:nVars
    pHF{ip} = plot(t_full,real(xr_hf_full(ip,:)),'r','LineWidth',1);
    hold on
    pLF{ip} = plot(t_full,real(xr_lf_full(ip,:)),'b','LineWidth',1);
end

for k = 1:nSplit %plot dotted lines between time splits
    t = lf_res{k}.t;
    subplot(2,1,1);
    plot([t(end) t(end)],get(gca, 'YLim'),'k:')
    subplot(2,1,2);
    plot([t(end) t(end)],get(gca, 'YLim'),'k:')
end

subplot(2,1,1)
xlim([t_full(1) t_full(end)]);
legend([pRaw{1},pRecon{1}],{'Raw Data','Full DMD Recon.'},'Location','best');
subplot(2,1,2)
xlim([t_full(1) t_full(end)]);
legend([pHF{1},pLF{1}],{'HF Recon.','LF Recon'},'Location','best')


%% Cluster HF Modes
rank_hf = 2; %rank of hf dynamics
allModes = [];
goodSplits = zeros(nSplit,1); %track which windows successfully yielded 2 hf modes
for k = 1:nSplit
    w = hf_res{k}.w;
    if (size(w,2) == rank_hf) & (isnan(w) == 0)
        allModes = [allModes w];
        goodSplits(k) = 1;
    end
end
allModes = [real(allModes); imag(allModes)].'; %separate real and imag for gmm
% hf_mode_gmm = fitgmdist(allModes,rank_hf);

[idx,~,~,clustDists] = kmeans(allModes,rank_hf);
clustLabels = zeros(rank_hf,nSplit);
splitCount = 0;
for k = 1:nSplit %sort modes in each window so all windows match up
    if goodSplits(k) == 0
        continue
    end
    splitCount = splitCount + 1; %just count non-error windows
    windDists = clustDists(rank_hf*(splitCount-1)+1:rank_hf*splitCount,:); %each row is that mode's distance to each of the centroids
    [~,naiveClass] = min(windDists,[],2);
    if length(unique(naiveClass)) == length(naiveClass)
        clustLabels(:,k) = naiveClass;
        continue %if all modes have been assigned to different clusters, we're done
    end
    
    windDists_pos = zeros(rank_hf,1);
    windDists_neg = zeros(rank_hf,1);
    
    for j = 1:rank_hf
        windDists_pos(j) = windDists(j,naiveClass(j)); %distances to positively-IDed modes
        windDists_neg(j) = sqrt(sum(windDists(j,:).^2) - windDists_pos(j).^2); %aggregate distances to negatively-IDed modes
    end
    
    windDists_ratio = windDists_pos./windDists_neg; %lower ratio = better confidence
    class_ratio = [naiveClass windDists_ratio];
    
    nDup = length(naiveClass) - length(unique(naiveClass));
    iter = 0;
    while nDup > 0
        iter = iter + 1;
        for j1 = 1:rank_hf
            for j2 = 1:rank_hf
                if ((naiveClass(j1) == naiveClass(j2)) && j1 ~= j2)
                    if class_ratio(j1) < class_ratio(j2) %kick the lower-confidence classification to a different label
                        if naiveClass(j2) == rank_hf 
                            naiveClass(j2) = 1;
                        else
                            naiveClass(j2) = naiveClass(j2) + 1;
                        end
                    else
                        if naiveClass(j1) == rank_hf 
                            naiveClass(j1) = 1;
                        else
                            naiveClass(j1) = naiveClass(j1) + 1;
                        end
                    end
                end
            end
        end
        nDup = length(naiveClass) - length(unique(naiveClass));
        if iter > 100
            disp('Infinite Loop')
            break
        end
    end
    clustLabels(:,k) = naiveClass;
end

sorted_modes_hf = cell(rank_hf,nSplit);
avg_modes_hf = zeros(nVars,rank_hf);
for k = 1:nSplit
    if goodSplits(k) == 0
        continue
    end
    w = hf_res{k}.w;
    b = hf_res{k}.b;
    Omega = hf_res{k}.Omega;
    t = hf_res{k}.t;
    
    if isempty(w)
        continue
    end
    
    w_sorted = w(:,clustLabels(:,k));
    
    %compute time-series projections onto hf modes
    bt = diag(b)*exp(Omega*t);
    hf_res{k}.bt = bt;
    
    for j = 1:rank_hf
        sorted_modes_hf{j,k} = w_sorted(:,j);
        avg_modes_hf(:,j) = avg_modes_hf(:,j) + w_sorted(:,j);
    end
    disp(k)
    disp(avg_modes_hf)
end
avg_modes_hf = avg_modes_hf/nnz(goodSplits);

for j = 1:rank_hf %normalize
    avg_modes_hf(:,j) = avg_modes_hf(:,j)/norm(avg_modes_hf(:,j));
end

x_proj_hf = avg_modes_hf.' * x_full;

figure
subplot(3,1,1)
p_hf_raw = plot(t_full,x_full,'LineWidth',1);
xlim([t_full(1),t_full(end)]);
title('Measurement Data')
legend('x_1','x_2','x_3','x_4','Location','eastoutside');
subplot(3,1,2)
p_hf_recon_re = plot(t_full,real(x_proj_hf),'LineWidth',1);
xlim([t_full(1),t_full(end)]);
title('Projection onto Avg. HF Modes (Real)')
legend('Re[b_1(t)]','Re[b_2(t)]','Location','eastoutside');
subplot(3,1,3)
p_hf_recon_im = plot(t_full,imag(x_proj_hf),'LineWidth',1);
xlim([t_full(1),t_full(end)]);
title('Projection onto Avg. HF Modes (Imag.)')
legend('Im[b_1(t)]','Im[b_2(t)]','Location','eastoutside');

%% Compute (Windowed) Projections onto HF Modes
rank_hf = 2;
window_size = size(hf_res{1}.x,2);
nGoodWind = nnz(goodSplits);
% b_hf_comb = zeros(1,nGoodWind*window_size);
% b_hf_comb_dt = zeros(1,nGoodWind*window_size);
b_hf_comb = zeros(1,nSteps);
b_hf_comb_dt = zeros(1,nSteps);
% t_good = zeros(1,nGoodWind*window_size);
splitCount = 0;
for k = 1:nSplit
    if goodSplits(k) == 0
        b_hf_comb_dt((k-1)*window_size + 1 : k*window_size) = NaN(1,window_size);
        b_hf_comb((k-1)*window_size + 1 : k*window_size) = NaN(1,window_size);
        continue
    end
    splitCount = splitCount + 1;
    x = hf_res{k}.x;
    w = hf_res{k}.w;
    t = hf_res{k}.t;
    bt = hf_res{k}.bt;
    w_sorted = w(:,clustLabels(:,k));
    
    %assume the 2 HF modes are complex conjugates, combine into one real
    %vector
    w_sorted_comb = (w_sorted(:,1) + w_sorted(:,2))/2;
    w_sorted_comb = abs(w_sorted_comb)/norm(w_sorted_comb); % strip any residual imaginary part & normalize
%     disp(w_sorted_comb)
%     phase_rot = angle(w_sorted(1,1));
%     w_sorted = w_sorted .* exp(-sqrt(-1)*phase_rot);
    
%     b = w_sorted.' * x;
%     b_hf_full(:,(k-1)*window_size + 1 : k*window_size) = b;

%     b_comb = w_sorted_comb.' * x;
    b_comb = real(sum(bt)); %real() is just to strip off any errant imaginary residue
%     b_hf_comb((splitCount-1)*window_size + 1 : splitCount*window_size) = b_comb;
    b_hf_comb((k-1)*window_size + 1 : k*window_size) = b_comb;
    
    b_comb_dt = zeros(size(b_comb));
    b_comb_dt(2:end-1) = (b_comb(3:end) - b_comb(1:end-2))/(2*(t(2)-t(1)));
    b_comb_dt(1) = (b_comb(2) - b_comb(1))/(t(2)-t(1));
    b_comb_dt(end) = (b_comb(end) - b_comb(end-1))/(t(2)-t(1));

%     b_hf_comb_dt((splitCount-1)*window_size + 1 : splitCount*window_size) = b_comb_dt;
    b_hf_comb_dt((k-1)*window_size + 1 : k*window_size) = b_comb_dt;

    t_good((splitCount-1)*window_size + 1 : splitCount*window_size) = t;
    %     disp(w_sorted)
%     disp([abs(w_sorted(:,1)) abs(w_sorted(:,2))])
end
figure
subplot(2,1,1)
plot(t_full,b_hf_comb,'LineWidth',1);
hold on
xlim([t_full(1),t_full(end)]);
title('Combined b_{HF}(t)')
hold on
subplot(2,1,2)
plot(t_full,b_hf_comb_dt,'LineWidth',1);
hold on
xlim([t_full(1),t_full(end)]);
title('\partial_t b_{HF}(t)')
hold on
% testom = 9;
% plot(t_full,3*abs(sin(testom*t_full)),'g')
% hold on
% plot(t_full,3*sin(testom*t_full).^2,'k')

% for k = 1:nSplit %plot dotted lines between time splits
%     t = hf_res{k}.t;
%     subplot(2,1,1);
%     plot([t(end) t(end)],get(gca, 'YLim'),'k','LineWidth',0.5)
%     subplot(2,1,2);
%     plot([t(end) t(end)],get(gca, 'YLim'),'k','LineWidth',0.5)
% end

hf_bt = [b_hf_comb; b_hf_comb_dt];
save('hf_sindy_data_2.mat','hf_bt','t_good','t_full','goodSplits','w_sorted_comb')