% Gaussian smoothing in time domain
% data is channels x time x trials (also works if just a single vector with time series)
function data = smooth_timeSeries(data,temp_smoothing,type)

    if any(~isnan(data(:)))

        if nargin<3
            type = 'Gauss';
        end

        if temp_smoothing>1

            if strcmp(type,'Gauss')
                w = gausswin(temp_smoothing);
                w = w./sum(w); % normalise
            elseif strcmp(type,'Uniform')
                w = (1/temp_smoothing)*ones(1,temp_smoothing);
            end

            if size(data,3)==1 && (size(data,1)==1 || size(data,2)==1)
                
                if strcmp(type,'Average')
                    data = smooth(data,temp_smoothing);
                elseif strcmp(type,'Gauss') || strcmp(type,'Uniform')
                    data = filter(w,1,data);
                end
                
            else
            
                for idx_channel=1:size(data,1)

                    for idx_trial=1:size(data,3)

                        if strcmp(type,'Average')
                            data(idx_channel,:,idx_trial) = smooth(data(idx_channel,:,idx_trial),temp_smoothing);
                        elseif strcmp(type,'Gauss') || strcmp(type,'Uniform')
                            data(idx_channel,:,idx_trial) = filter(w,1,data(idx_channel,:,idx_trial));
                        end

                    end

                end  
                
            end

        end  

    end
    


end