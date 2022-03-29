% baseline correct time series for a specified interval

function data = baseline_correct_timeSeries(data,interval)


[n_chan, ~, n_trials] = size(data);

for idx_chan=1:n_chan
    
    for idx_trial=1:n_trials
        
        baseline_mean = nanmean(data(idx_chan,interval(1):interval(2),idx_trial));
        
        data(idx_chan,:,idx_trial) = data(idx_chan,:,idx_trial) - baseline_mean;
        
    end
    
end


end
