% function to make RSA for time series

function [RSM, beta_all] = mk_RSA_timeseries(labStm,n_time,data,X,idx_chan_use,plot_chan)
    
    n_states = size(X,2);

    RSM = nan(size(X,2)*(size(X,2)-1)/2,n_time);
    
    beta_all = nan(size(X,2),size(data,1),n_time);
    
    X = [X ones(size(X,1),1)];
    
    for idx_time=1:n_time
        
        Y = squeeze(data(:,idx_time,:));
        
        cm = eye(size(X,2));
        
        beta = ols_PS(Y',X,cm);
        
        if size(beta,1)>n_states
            beta_all(:,:,idx_time) = beta(1:end-1,:);
        else
            beta_all(:,:,idx_time) = beta;
        end
        
        % pre-whiten (based on rsa-toolbox code!):        
        res  = Y' - X*beta;
        
        Sw_reg = covdiag_PS(res,size(res,1)-1,'shrinkage');
        
        [V,L] = eig(Sw_reg);       % Faster and numerically more stable than Sw_hat.^-1/2
        l     = diag(L);
        sq    = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        u_hat = beta*sq;
        
        if size(beta,1)>n_states
            RSM_trial =  corr(u_hat(1:end-1,:)');
        else
            RSM_trial =  corr(u_hat');
        end
        
        RSM(:,idx_time) = RSM_trial(tril(true(length(RSM_trial)),-1));
        
    end
    
    if plot_chan
       beta_all_plot = nanmean(beta_all(:,:,71:150),3); % mean over given time interval
       beta_all_plot = squeeze(nanmean(beta_all_plot(:,:),1)); % mean over stims
       
       zPlotSens_PS( beta_all_plot, idx_chan_use)
    end
    
end