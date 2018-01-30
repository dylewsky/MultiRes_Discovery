function [res_list, kmList, mr_res] = example2_tdDMD(x,TimeSpan,wSteps,nSplit,r,use_last_freq,thresh_pct,nComponents)
    % addpath('optdmd-master');
    %addpath(genpath(fullfile('..','..','Optimized DMD','optdmd-master')));
    % load('../raw_data_2_hiRes.mat');


    % delaySteps = 1000;
    % nDelay = 5;
    

    % r = size(x,1); %rank to fit w/ optdmd
    % r = 12;
    imode = 1;
    %  imode = 1, fit full data, slower
    %  imode = 2, fit data projected onto first r POD modes
    %      or columns of varargin{2} (should be at least r
    %      columns in varargin{2})
    nVars = size(x,1);

    
    nSteps = wSteps*nSplit;

    % use_last_freq = 0; %use previous spectrum as guess for next window of dmd


    primeList = 1;

    %% OptDMD

    mr_res = cell(length(primeList),nSplit);
    res_list = [];
    for pn = 1:length(primeList)
        xSample = x(:,1:nSteps);
        tSample = TimeSpan(1:nSteps);

        nHold = 0;
        xL = reshape(xSample, nVars, wSteps, nSplit);
        tL = reshape(tSample, 1, wSteps, nSplit);

        res_list = [res_list; pn, nSplit, wSteps];
        clear('e_init')
        for k = 1:nSplit
            xT = xL(:,:,k);
            tT = tL(:,:,k);
            mr_res{pn,k}.x = xT;
            mr_res{pn,k}.t = tT;
            try
                c = mean(xT,2);
                xT = xT - repmat(c,1,size(xT,2));
                t_start = tT(1);
                tT = tT - t_start;
                if (exist('e_init','var')) && (use_last_freq == 1)
                    [w, e, b] = optdmd(xT,tT,r,imode,[],e_init);
                else
                    [w, e, b] = optdmd(xT,tT,r,imode);
                end
            catch ME
                disp(['Error on pn = ' num2str(pn) ', k = ' num2str(k)]);
                continue
            end
    %             Omega = log(diag(e))/(TimeSpan(2) - TimeSpan(1));
            e_init = e;
            mr_res{pn,k}.c = c;
            mr_res{pn,k}.t_start = t_start;
            mr_res{pn,k}.w = w;
            mr_res{pn,k}.Omega = e;
            mr_res{pn,k}.b = b;
    %             mr_res{pn,k}.Omega = Omega;
        end
    end


    %% Cluster Frequencies
    if exist('mr_res','var') == 0
        load('mr_res_2_td.mat');
        load('res_list_2_td.mat');
    end

    % gmmList = cell(size(res_list,1),1);
    kmList = zeros(size(res_list,1),nComponents);

    for q = 1:size(res_list,1)
    % for q = 16
        pn = primeList(q);
        om_spec = zeros(nVars,nSteps);

        all_om = [];
        all_om_grouped = [];
        for k = 1:nSplit

            try
                Omega = mr_res{pn,k}.Omega;
            catch ME
                disp(['Error on q = ' num2str(q) ', k = ' num2str(k)]);
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

        log_om_sq = log(all_om_sq);

        [~, C] = kmeans(all_om_sq,nComponents,'Distance','cityblock','Replicates',5);

        [C,sortInd] = sort(C);

        kmList(q,:) = C.';

        for k = 1:nSplit
            try
                omega = mr_res{pn,k}.Omega;
            catch ME
                continue
            end
            om_sq = omega.*conj(omega);

            om_sq_dists  = (repmat(C.',r,1) - repmat(om_sq,1,nComponents)).^2;
            [~,om_class] = min(om_sq_dists,[],2);


            mr_res{pn,k}.om_class = om_class;

        end

    end

%     figure
%     mColors = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5; 0.33 0.33 0.33; 0 0 0; 0.66 0.66 0.66; 1 1 0; 1 0 1; 0 1 1];
%     for k = 1:size(all_om_grouped,1)
%         for jk = 1:size(all_om_grouped,2)
%             plot(real(all_om_grouped(k,jk)), imag(all_om_grouped(k,jk)),'.','Color',mColors(jk,:),'MarkerSize',5)
%             hold on
%         end
%     end
%     xlim([-40 40]);
%     ylim([-80 80]);
%     save('mr_res_2_td.mat', 'mr_res');
%     save('res_list_2_td.mat', 'res_list');
%     % save('gmm_list_2_td.mat', 'gmmList');
%     save('km_list_2_td.mat','kmList');
   
end
