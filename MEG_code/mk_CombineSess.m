% function to compute nanmean over sessions to summarise into before and
% after rest
function seq_data_sess = mk_CombineSess(seq_data,n_Sess_Combine)

    if any(~isnan(seq_data(:)))

        if nargin<2
            n_Sess_Combine = 3;
        end

        n_sub   = size(seq_data,1);
        n_time  = size(seq_data,2);
        n_sess  = size(seq_data,3);
        n_shuff = size(seq_data,4);
        n_from   = size(seq_data,5);
        n_to     = size(seq_data,6);

        n_sess_new = n_sess/n_Sess_Combine; % number of trials per session

        seq_data_sess = nan(n_sub,n_time,n_sess_new,n_shuff,n_from,n_to);

        sess_bounds = 1:n_Sess_Combine:n_sess;

        for idx_sess=1:n_sess_new
            seq_data_sess(:,:,idx_sess,:,:,:) = nanmean(seq_data(:,:,sess_bounds(idx_sess):sess_bounds(idx_sess)+(n_Sess_Combine-1),:,:,:),3);
        end
        
    else
        
        seq_data_sess = seq_data;
        
    end

end