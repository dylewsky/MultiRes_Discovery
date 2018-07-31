clear variables; close all; clc

load('mwDMD_mr_res_i2.mat');
load('mwDMD_params.mat');
load('mwDMD_sep_recon.mat');

forecast_length = 5*wSteps;
rng(111); %seed
test_windows = randperm(nSlide-1000) + 500;
test_windows = sort(test_windows(1:3));

suppress_growth = 1; %eliminate all positive real components to frequencies

poly_fit_order = 5;

dt = mr_res{1}.t(2) - mr_res{1}.t(1);
%% DMD vs. Polynomial vs HAVOK Forecast
for j = 1:length(test_windows)
% for j = 3
    t_test = mr_res{test_windows(j)}.t_start + dt*(0:forecast_length-1);
    steps_test = round(t_test/dt);
    
    w = mr_res{test_windows(j)}.w;
    b = mr_res{test_windows(j)}.b;
    
    Omega = mr_res{test_windows(j)}.Omega;
    if suppress_growth == 1
        Omega(real(Omega) > 0) = sqrt(-1) * imag(Omega(real(Omega) > 0));
    end
    om_class = mr_res{test_windows(j)}.om_class;
    t = mr_res{test_windows(j)}.t;
    c = mr_res{test_windows(j)}.c;
    t_start = mr_res{test_windows(j)}.t_start;
    tShift = t-t(1);
    
    figure
    for k = 1:nComponents
        
        % DMD Forecast
        
        xr_comp = w(:, om_class == k)*diag(b(om_class == k))*exp(Omega(om_class == k)*(t_test-t_start));
        if k == 1 % constant shift gets put in LF recon
            xr_comp = xr_comp + c; 
        end
        xr_comp = real(xr_comp);
        
        subplot(nComponents,4,4*(k-1)+1)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['DMD Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(t_test,xr_comp,'r-','DisplayName','1-Window DMD Recon')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
%         legend('Location','eastoutside');

        % Polynomial Forecast
        
        t_train = t_test(1:wSteps);
        xr_train = xr_sep{k}(:,steps_test(1:wSteps));
        xr_mean = mean(xr_train,2);
        xr_train = xr_train - xr_mean;
        pc = cell(nVars,1); %polynomial fit coefficients for each variable
        poly_forecast = cell(nVars,1);
        for vi = 1:nVars
            pc{vi} = polyfit(t_train,xr_train(vi,:),poly_fit_order);
            poly_forecast{vi} = polyval(pc{vi},t_test) + xr_mean(vi);
        end
        yBounds = ylim;
        subplot(nComponents,4,4*(k-1)+2)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Polynomial Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        for vi = 1:nVars
            plot(t_test,poly_forecast{vi},'r-','DisplayName','1-Window Poly Recon')
            hold on
        end
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
        
        %Global HAVOK Forecast
        
        load(['MR_HAVOK_model_' num2str(k) '.mat']);
        x0 = zeros(nDelay*size(xr_sep{k},1),1);
%         t0 = zeros(nDelay,1);
        for nd = 1:nDelay
            delayInds = ((nd-1)*size(xr_sep{k},1) + 1) : nd*size(xr_sep{k},1);
            x0(delayInds) = xr_sep{k}(:,steps_test(wSteps) - (nd-1)*delaySteps);
%             t0(nd) = tspan(steps_test(wSteps) - (nd-1)*delaySteps);
        end
        v0 = S^(-1) * U.' * x0; %transform into eigenbasis used by HAVOK model
        v0 = v0(1:r-1); %...up to rank used by HAVOK model

        [tHAV, vHAV] = ode45(@(tHAV,vHAV)HAVOK_rhs(vHAV,A),t_test(wSteps-1:end),v0);
        xHAV = U(:,1:r-1) * S(1:r-1,1:r-1) * vHAV.';
        xHAV = xHAV(1:nVars,:);
        subplot(nComponents,4,4*(k-1)+3)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Global HAVOK Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(tHAV,xHAV,'r-','DisplayName','HAVOK Forecast')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
        
        % Local HAVOK Forecast
        
        % EIGEN-TIME DELAY COORDINATES
        Hloc = zeros(nVars*nDelay,length(steps_test)-(nDelay-1)*delaySteps);
        for hk=1:nDelay
            delayInds = ((hk-1)*nVars + 1) : hk*nVars;
            Hloc(delayInds,:) = xr_sep{k}(:,(hk-1)*delaySteps+1:(length(steps_test)-(nDelay)*delaySteps + hk*delaySteps));
        end

        [Uloc,Sloc,Vloc] = svd(Hloc,'econ'); % Eigen delay coordinates

        % COMPUTE DERIVATIVES (4TH ORDER CENTRAL DIFFERENCE)
        dVloc = zeros(length(Vloc)-5,r);
        for i=3:length(Vloc)-3
            for hk=1:r
                dVloc(i-2,hk) = (1/(12 * dt)) * (-Vloc(i+2,hk)+8 * Vloc(i+1,hk)-8 * Vloc(i-1,hk)+Vloc(i-2,hk));
            end
        end
        % trim first and last two that are lost in derivative
        Vloc = Vloc(3:end-3,1:r);

        % BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
        XiLoc = Vloc\dVloc;
        Aloc = XiLoc(1:r-1,1:r-1)';
        Bloc = XiLoc(end,1:r-1)';
        
        v0loc = Sloc^(-1) * Uloc.' * x0; %transform into eigenbasis used by HAVOK model
        v0loc = v0loc(1:r-1); %...up to rank used by HAVOK model

        [tHAVloc, vHAVloc] = ode45(@(tHAVloc,vHAVloc)HAVOK_rhs(vHAVloc,Aloc),t_test(wSteps-1:end),v0loc);
        xHAVloc = Uloc(:,1:r-1) * Sloc(1:r-1,1:r-1) * vHAVloc.';
        xHAVloc = xHAVloc(1:nVars,:);
        subplot(nComponents,4,4*(k-1)+4)
        plot(t_test,xr_sep{k}(:,steps_test),'k-','DisplayName','Full DMD Recon.')
        title(['Local HAVOK Forecast at t = ' num2str(t_test(1)) ': Component ' num2str(k)]);
        hold on
        plot(tHAVloc,xHAVloc,'r-','DisplayName','HAVOK Forecast')
        hold on
        plot([t_test(wSteps) t_test(wSteps)],ylim,'k--','DisplayName','Boundary Between Window and Future')
        hold off
        xlim([t_test(1) t_test(end)])
        ylim(yBounds)
    end
end
