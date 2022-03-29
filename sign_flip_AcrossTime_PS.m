% This function does sign-flip permutation, which gives a prob distribution
% for a DIFFERENCE value having a mean of zero
% 
% INPUT:
% data: a data vector of mean differences
% n_perm: how many permutations
% 
% OUTPUT
% up_bound: upper bound of 95% of permutation distribution
% down_bound: lower bound of 95% of permutation distribution
% emp_percentile: percentile of empirical (observed) difference

function [up_bound, down_bound, emp_percentile, sig_t_maxClus, diff_values_output] = sign_flip_AcrossTime_PS(data,n_perm,diff_values,perc_bounds,do_tvals,P_t,one_sided)

if nargin<2
    n_perm = 1000;
%     perc_bounds = [97.5,2.5];
    perc_bounds = [100,0];
    do_tvals    = false;
    P_t         = 0.05;
    one_sided   = false;
elseif nargin<4
%     perc_bounds = [97.5,2.5];
    perc_bounds = [100,0];
    do_tvals    = false;
    P_t         = 0.05;
    one_sided   = false;
elseif nargin<5
    do_tvals    = false;
    P_t         = 0.05;
    one_sided   = false;
elseif nargin==0
    error('need to provide some data!')
end

if do_tvals
    t_values = nan(size(data,1),size(data,2),n_perm);
    if ~one_sided
        crit_t   = tinv(1-P_t/2,size(data,1)-1);
    else
        crit_t   = tinv(1-P_t,size(data,1)-1);
    end
    sig_t_maxClus  = nan(1,n_perm);
else
    sig_t_maxClus  = [];
end

if nargin<3 || isempty(diff_values)
    
    diff_values = nan(size(data,1),size(data,2),n_perm);
    
    for idx_perm=1:n_perm
        if idx_perm==1
            diff_values(:,:,idx_perm) = data;
        else
            diff_values(:,:,idx_perm) = data.*repmat(randsample([-1 1],size(data,1),true)',1,size(data,2));
        end
        
        if do_tvals
            
            t_vals_perm = mean(diff_values(:,:,idx_perm))./(std(diff_values(:,:,idx_perm))/sqrt(size(data,1)));
            
            if strcmp(one_sided,'one_sided_neg')
                exceed_t_pos = 0;
            else
                exceed_t_pos = t_vals_perm>crit_t;
                exceed_t_pos = max(accumarray(nonzeros((cumsum(~exceed_t_pos)+1).*exceed_t_pos),1));
                if isempty(exceed_t_pos)
                    exceed_t_pos = 0;
                end
            end
            
            if strcmp(one_sided,'one_sided_pos')
                exceed_t_neg = 0;
            else
                exceed_t_neg = t_vals_perm<-crit_t;
                exceed_t_neg = max(accumarray(nonzeros((cumsum(~exceed_t_neg)+1).*exceed_t_neg),1));
                if isempty(exceed_t_neg)
                    exceed_t_neg = 0;
                end
            end
            
            sig_t_maxClus(idx_perm) = max([exceed_t_pos exceed_t_neg]);
            
        end
        

    end
    
end

diff_values_output = diff_values;

diff_values = squeeze(nanmean(diff_values,1));

up_bound   = prctile(diff_values(:,2:end),perc_bounds(1),2); % exclude true data
down_bound = prctile(diff_values(:,2:end),perc_bounds(2),2); % exclude true data

emp_percentile = nan(size(diff_values,1),1);
for idx_lag = 1:size(diff_values,1)
    emp_percentile(idx_lag) = (sum(diff_values(idx_lag,:) < diff_values(idx_lag,1)) + 0.5*sum(diff_values(idx_lag,:) == diff_values(idx_lag,1))) / length(diff_values(idx_lag,:)); 
end


end