% function to compute nanmean over trials to summarise the six sessions
function seq_data_sess = mk_MeanSess(seq_data,n_trials_perSess)

if nargin<2
    n_trials_perSess = 48;
end

if any(~isnan(seq_data(:)))

    n_sub    = size(seq_data,1);
    n_time   = size(seq_data,2);
    n_trials = size(seq_data,3);
    n_shuff  = size(seq_data,4);
    n_from   = size(seq_data,5);
    n_to     = size(seq_data,6);    

    n_sess = n_trials/n_trials_perSess; % number of trials per session

    seq_data_sess = nan(n_sub,n_time,n_sess,n_shuff,n_from,n_to);

    sess_bounds = 1:n_trials_perSess:n_trials;

    for idx_sess=1:n_sess
        seq_data_sess(:,:,idx_sess,:,:,:) = nanmean(seq_data(:,:,sess_bounds(idx_sess):sess_bounds(idx_sess)+(n_trials_perSess-1),:,:,:),3);
    end

else
    
    seq_data_sess = seq_data;
    
end

end