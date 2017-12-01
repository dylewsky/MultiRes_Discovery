clear variables;
load('res_list.mat');
load('mr_res.mat');
load('raw_data.mat');


sepTrial = 16; %index of trial with window size that successfully separated scales
lfScale = 1;
hfScale = 2; % 1 = low-frequency, 2 = high-frequency
% sepModes = 2; % # of modes 

nVars = size(x,1);
nSteps = 2^14;

j = res_list(sepTrial,2);
pn = res_list(sepTrial,1);
nSplit = 2^(j-1);

hf_res = cell(nSplit, 1);
lf_res = cell(nSplit, 1);
for k = 1:nSplit
    om_class = mr_res{pn,j,k}.om_class;
    hf_res{k}.t = mr_res{pn,j,k}.t;
    hf_res{k}.x = mr_res{pn,j,k}.x(om_class == hfScale,:);
    hf_res{k}.w = mr_res{pn,j,k}.w(:,om_class == hfScale);
    hf_res{k}.Omega = mr_res{pn,j,k}.Omega(om_class == hfScale);
    hf_res{k}.b = mr_res{pn,j,k}.b(om_class == hfScale);
    hf_res{k}.om_post = mr_res{pn,j,k}.om_post(om_class == hfScale,:);
    
    lf_res{k}.t = mr_res{pn,j,k}.t;
    lf_res{k}.x = mr_res{pn,j,k}.x(om_class == lfScale,:);
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
    xr_window = w*diag(b)*exp(Omega*t);
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
    pRaw{ip} = plot(t_full,x_full(ip,:),'k','LineWidth',2);
    hold on
    pRecon{ip} = plot(t_full,real(xr_full(ip,:)),'g','LineWidth',2);
end
subplot(2,1,2)
pHF = cell(nVars,1);
pLF = cell(nVars,1);
for ip = 1:nVars
    pHF{ip} = plot(t_full,real(xr_hf_full(ip,:)),'r','LineWidth',2);
    hold on
    pLF{ip} = plot(t_full,real(xr_lf_full(ip,:)),'b','LineWidth',2);
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
for k = 1:nSplit
    w = hf_res{k}.w;
    allModes = [allModes w];
end
allModes = [real(allModes); imag(allModes)].'; %separate real and imag for gmm
% hf_mode_gmm = fitgmdist(allModes,rank_hf);

[idx,~,~,clustDists] = kmeans(allModes,rank_hf);
clustLabels = zeros(rank_hf,nSplit);
for k = 1:nSplit
    windDists = clustDists(rank_hf*(k-1)+1:rank_hf*k,:); %each row is that mode's distance to each of the centroids
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
    w = hf_res{k}.w;
    
    w_sorted = w(:,clustLabels(:,k));

    for j = 1:rank_hf
        sorted_modes_hf{j,k} = w_sorted(:,j);
        avg_modes_hf(:,j) = avg_modes_hf(:,j) + w_sorted(:,j);
    end
end
avg_modes_hf = avg_modes_hf/nSplit;

for j = 1:rank_hf %normalize
    avg_modes_hf(:,j) = avg_modes_hf(:,j)/norm(avg_modes_hf(:,j));
end

x_proj_hf = avg_modes_hf.' * x_full;

figure
plot(t_full,x_full,'k:',t_full,x_proj_hf,'r-','LineWidth',2)




