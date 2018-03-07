function [mr_res, km_centroids] = run_optDMD(x,TimeSpan,r,wSteps,nSplit,stepSize,nComponents,thresh_pct,use_last_freq,use_median_freqs)
    imode = 1; %parameter for optdmd code

    nSteps = wSteps * nSplit;
    nVars = size(x,1);

    nSlide = floor((nSteps-wSteps)/stepSize);
    
    global_SVD = 1; %use global SVD modes for each DMD rather than individual SVD on each window

    if global_SVD == 1
        [u,~,v] = svd(x,'econ');
        u = u(:,1:r); %just first r modes
        v = v(:,1:r);
    end

    corner_sharpness = 16; %higher = sharper corners
    lv_kern = tanh(corner_sharpness*(1:wSteps)/wSteps) - tanh(corner_sharpness*((1:wSteps)-wSteps)/wSteps) - 1;

    mr_res = cell(nSlide,1);
    for k = 1:nSlide
        sampleStart = stepSize*(k-1) + 1;
        sampleSteps = sampleStart : sampleStart + wSteps - 1;
        xSample = x(:,sampleSteps);
        tSample = TimeSpan(sampleSteps);

        mr_res{k}.x = xSample;
        mr_res{k}.t = tSample;

        c = mean(xSample,2); %subtract off mean before rounding corners
        xSample = xSample - repmat(c,1,size(xSample,2));

        xSample = xSample.*repmat(lv_kern,nVars,1); %round off corners

        t_start = tSample(1);
        tSample = tSample - t_start;
        if global_SVD == 0
            [u,~,~] = svd(xSample,'econ');
            u = u(:,1:r);
        end
        if (exist('e_init','var')) && (initialize_artificially == 1)
            [w, e, b] = optdmd(xSample,tSample,r,imode,[],e_init,u);
        else
            [w, e, b] = optdmd(xSample,tSample,r,imode,[],[],u);
        end
        if use_last_freq == 1
            e_init = e;
        end
        if use_median_freqs == 1
            [eSq, eInd] = sort(e.*conj(e)); %match order to that of freq_meds
            freq_angs = angle(e(eInd));
            e_init = exp(sqrt(-1)*freq_angs) .* freq_meds;
        end
        mr_res{k}.w = w;
        mr_res{k}.Omega = e;
        mr_res{k}.b = b;
        mr_res{k}.c = c;
        mr_res{k}.t_start = t_start;
    end
    

    %% Cluster Frequencies
%     close all;
    if exist('mr_res','var') == 0
        if nonlinear == 0
            try
                load('mwDMD_mr_res_i2_linear.mat');
            catch ME
                load('mwDMD_mr_res_linear.mat');
            end
        else
            try
                load('mwDMD_mr_res_i2.mat');
            catch ME
                load('mwDMD_mr_res.mat');
            end
        end
    end

    % nBins = 64;



    all_om = [];
    all_om_grouped = [];
    for k = 1:nSlide
        try
            Omega = mr_res{k}.Omega;
        catch ME
            disp(['Error on k = ' num2str(k)]);
            continue
        end
        all_om = [all_om; Omega];
        [~,imSortInd] = sort(imag(Omega));
        all_om_grouped = [all_om_grouped; Omega(imSortInd).'];
    end

    all_om_sq = conj(all_om) .* all_om;
    all_om_sq = sort(all_om_sq);
    all_om_sq_thresh = all_om_sq(floor(thresh_pct*length(all_om_sq)));
    all_om_sq = all_om_sq(1:floor(thresh_pct*length(all_om_sq)));


    [~, km_centroids] = kmeans(all_om_sq,nComponents,'Distance','cityblock','Replicates',5);

    [km_centroids,sortInd] = sort(km_centroids);

    for k = 1:nSlide
        try
            omega = mr_res{k}.Omega;
        catch ME
            continue
        end
        om_sq = omega.*conj(omega);

        om_sq_dists  = (repmat(km_centroids.',r,1) - repmat(om_sq,1,nComponents)).^2;
        [~,om_class] = min(om_sq_dists,[],2);

        mr_res{k}.om_class = om_class;
    end
    if exist('mr_res','var') == 0
        if nonlinear == 0
            if use_median_freqs == 0
                save('mwDMD_mr_res_linear.mat', 'mr_res');
            else
                save('mwDMD_mr_res_i2_linear.mat', 'mr_res');
            end
        else
            if use_median_freqs == 0
                save('mwDMD_mr_res.mat', 'mr_res');
            else
                save('mwDMD_mr_res_i2.mat', 'mr_res');
            end
        end
    end
end