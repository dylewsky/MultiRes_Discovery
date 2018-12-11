load('recursive_params.mat')
load('Three_Body_Data_Cartesian.mat')
nVars = size(pos,1);
tspan_orig = tspan;
tot_components = sum(rec_nComponents) - sum(cell2mat(rec_composite_components));

all_tBounds = zeros(tot_components,2);
all_dt = zeros(tot_components,1);

figure
subplot(tot_components+1,1,1)
c_ind = 1;
plot(tspan,pos,'k')
ylabel('Full Signal')
for nr = 1:n_recursions
    load(['mwDMD_sep_recon_recursion_' num2str(nr,'%02.f') '.mat']);
    nComponents = rec_nComponents(nr);
    all_tBounds(nr,:) = [tspan(1) tspan(end)];
    all_dt(nr) = tspan(2)-tspan(1);
    for k = 1:nComponents
        if ismember(k,rec_composite_components{nr})
            continue
        end
        c_ind = c_ind + 1;
        subplot(tot_components+1,1,c_ind)
        plot(tspan,xr_sep{k})
        ylabel(['Pass ' num2str(nr) ' Component ' num2str(k)])
    end
end


%% 
tot_tBounds = [max(all_tBounds(:,1)) min(all_tBounds(:,2))];
min_dt = min(all_dt);
tspan_comb = tot_tBounds(1):min_dt:tot_tBounds(end);

sum_recon = zeros(nVars,length(tspan_comb));

c_ind = 1;
for nr = 1:n_recursions
    nComponents = rec_nComponents(nr);
    load(['mwDMD_sep_recon_recursion_' num2str(nr,'%02.f') '.mat']);
    for k = 1:nComponents
        if ismember(k,rec_composite_components{nr})
            continue
        end
        c_ind = c_ind + 1;
        for q = 1:nVars
            sum_recon(q,:) = sum_recon(q,:) + interp1(tspan, ...
                xr_sep{k}(q,:),tspan_comb);
        end
    end
end

if all_dt(1) ~= min_dt
    disp('Time Step Mismatch!')
end

pos_trunc = pos(:,1:length(tspan_comb));

scale_loss = sum_recon./pos_trunc;
scale_loss = mean(scale_loss);
scale_loss_median = median(scale_loss);
% figure
% plot(tspan_comb,scale_loss);
% ylim([-0.1 1.5])

figure
plot(tspan_orig,pos,'k','LineWidth',1.5)
hold on
plot(tspan_comb,(1/scale_loss_median)*sum_recon,'b','LineWidth',1)
