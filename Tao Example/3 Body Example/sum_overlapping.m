function [t_join,filt_join] = sum_overlapping(all_t,all_x)
% Sums a group of time series (t/x vector pairs) which partially overlap
% with one another
% all_t and all_x are matrices of the t and x vectors (concatenated in
% either direction, but should have same dimensions)
    all_t = reshape(all_t,numel(all_t),1);
    all_x = reshape(all_x,numel(all_x),1);
    
    [G,t_join] = findgroups(all_t);
    filt_join = splitapply(@sum,all_x,G);
end

