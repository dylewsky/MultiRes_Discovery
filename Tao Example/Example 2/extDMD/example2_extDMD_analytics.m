load('mr_res_2_ext.mat');

qBest = 5;
nSplit = 128;
r = length(mr_res{qBest,1}.Omega);
nFeat = size(mr_res{qBest,1}.x,1);

allModes = zeros(nFeat,r,nSplit);

for k = 1:nSplit
    try
        omega = mr_res{qBest,k}.Omega;
    catch ME
        disp(['Rank Deficient at k = ' num2str(k)])
        continue
    end
    om_sq = omega.*conj(omega);
    [om_sq, sind] = sort(om_sq);
    w = mr_res{qBest,k}.w;
    if nnz(w == Inf) ~= 0
        disp('Infinite-Valued Mode');
        continue
    elseif nnz(isnan(w)) ~= 0
        disp('NaN in Mode');
        continue
    end
    w = w(:,sind);
    allModes(:,:,k) = w;
end

absModes = abs(allModes);
meanModes = nansum(absModes,3);