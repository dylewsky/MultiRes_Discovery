function mr_res_sds = get_mr_res_stats(mr_res,nc)
    nSlide = length(mr_res);
    all_omega = [];
    all_omega_ang = [];
    all_b = [];
    all_w = [];
    
    for k = 1:nSlide
        if nnz(mr_res{k}.om_class == nc) ~= 2
            continue;
        end
        all_omega = [all_omega; mr_res{k}.Omega(mr_res{k}.om_class == nc)];
        all_omega_ang = [all_omega_ang; ...
            [abs(mr_res{k}.Omega(mr_res{k}.om_class == nc))...
            abs(angle(mr_res{k}.Omega(mr_res{k}.om_class == nc)))]];
        all_b = [all_b; mr_res{k}.b(mr_res{k}.om_class == nc)];
        this_w = mr_res{k}.w(:,mr_res{k}.om_class == nc);
        all_w = [all_w this_w(:,1)]; %just keep one of conjugate pair
    end
    all_w_ang = [abs(all_w) angle(all_w)];
    
    mr_res_sds.omega =  [std(real(all_omega)) std(abs(imag(all_omega)))];
    mr_res_sds.omega_ang = [std(all_omega_ang(:,1)) std(all_omega_ang(:,2))];
    mr_res_sds.b = [std(real(all_b)) std(abs(imag(all_b)))];
    mr_res_sds.w_ang = [std(all_w_ang(:,1),0,2) std(abs(all_w_ang(:,2)),0,2)];
    
end

