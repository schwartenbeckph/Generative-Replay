
function [betas, intercepts, corr_betas, amount_null, n_nonZ_chann, betas_11, ...
          intercepts_11, idx_chan, betas_struct, intercepts_struct, lab_struct, ...
          preds_obt,labStm,preds_obt_struct,labStm_struct] = obtain_betas(scan_result_path,name_file,temp_smoothing,ts_best,include_null,optimise_null,L1_use,do_normalise,n_null,threshold_corr,plot_data,baseline_correct,Loc_11,which_chan,check_accuracy)

    preds_obt_struct = NaN;
    labStm_struct    = NaN;
      
    if ~Loc_11
        betas_11          = nan;
        intercepts_11     = nan;
        betas_struct      = nan;
        intercepts_struct = nan;
        lab_struct        = nan;
    end
    
    if nargin<15
        check_accuracy = false;
    end

    %%%%% obtain imaging data (localiser and task) and smooth data %%%%%
    load(fullfile(scan_result_path,name_file)) % this loads 'data' and 'stimlabel'
    % data localiser is 201, the first 50 is prior to onset

    % take specific channels?
    channel_names = char(channel_names');
    
    if strcmp(which_chan,'all')
        idx_chan = true(size(data,1),1);
    elseif strcmp(which_chan,'Occ')
        idx_chan = channel_names(:,3)=='O';
    elseif strcmp(which_chan,'Front')
        idx_chan = channel_names(:,3)=='F'; 
    elseif strcmp(which_chan,'Temp')
        idx_chan = channel_names(:,3)=='T'; 
    elseif strcmp(which_chan,'Par')
        idx_chan = channel_names(:,3)=='P'; 
    elseif strcmp(which_chan,'Cen')
        idx_chan = channel_names(:,3)=='C';         
    end
    
    data = data(idx_chan,:,:);
    
    %%%%% smooth data (optional) %%%%%
    data = smooth_timeSeries(data,temp_smoothing);
    
    if baseline_correct
       data = baseline_correct_timeSeries(data,[41,50]); 
%        data = baseline_correct_timeSeries(data,[1,50]); 
    end

    labStm    = zeros(1,length(stimlabel));
    allconds  = unique(stimlabel);
    num_label = 1:length(allconds);

    for idx_num_label=1:max(num_label)
        labStm(strcmp(stimlabel,allconds{idx_num_label})) = idx_num_label;
    end

    %%%%% train classifiers on training data %%%%%
    nulldata = squeeze(data(:,1,:))';
    data     = squeeze(data(:,ts_best,:))';       

    % normalise data
    idx_Nonan     = ~isnan(data(:,1));
    dataLI        = data(idx_Nonan,:);  
    labStm        = labStm(idx_Nonan);
    
    if do_normalise
        dataLI = mk_normalise(dataLI);
    end
    
    if include_null
        idx_Nonan_null = ~isnan(nulldata(:,1));
        nullDataLI     = nulldata(idx_Nonan_null,:);
        if do_normalise
            nullDataLI = mk_normalise(nullDataLI);
        end
    else
        nullDataLI = [];
    end
    
    [betas, intercepts, amount_null] = train_classifiers(dataLI,nullDataLI,labStm,L1_use,optimise_null,[],[],false);    
    
    if Loc_11
        [betas_11, intercepts_11, amount_null] = train_classifiers(dataLI,nullDataLI,labStm,L1_use,optimise_null,[],[],true);
    end
    
    if plot_data
        figure,imagesc([dataLI; nullDataLI]),title(sprintf('Training Data')),colorbar
        figure,imagesc(betas),title(sprintf('Sub: %d, L1: %.3f',1,L1_use)),colorbar
        fprintf('Number non-zero channels per classifier is %d %d %d %d.\n',sum(betas>0))
    end
    
    if Loc_11
        r = corr([squeeze(betas_11(:,1,[2 3 4])) squeeze(betas_11(:,2,[1 3 4])) ...
                  squeeze(betas_11(:,3,[1 2 4])) squeeze(betas_11(:,4,[1 2 3]))]);
              
        n_nonZ_chann = sum(betas_11>0);
    else
        r = corr(betas);
        
%         fprintf('Correlations are %f %f %f %f %f %f.\n',r(tril(true(length(r)),-1))')
    
        fprintf('Number non-zero channels per classifier is %d %d %d %d.\n',sum(betas>0))

        n_nonZ_chann = sum(betas>0);
    end    
    
    fprintf('Mean correlation is %f.\n',mk_meanCorr_PS(r(tril(true(length(r)),-1))'))
    
    if any(isnan(r(tril(true(length(r)),-1))))
        error('Wrong classifiers')
    end
    
    corr_betas = r(tril(true(length(r)),-1))';
    
    preds_obt  = 1 ./ (1 + exp(-(dataLI*betas + intercepts)));        
    [~,plabs]  = max(preds_obt,[],2);        
    accur_mean = mean(plabs==labStm');
    fprintf('The mean accuracy of classifiers is %f.\n',accur_mean) 
    
    if check_accuracy            
        
        plot_identifyabilty_PS(preds_obt,labStm',[0 0.6])
        
        if (exist('lab_struct','var') && any(any(~isnan(lab_struct)))) || strcmp(name_file,'Data_inference.mat')
            if length(unique(labStm_struct))<8
                plot_identifyabilty_PS(preds_obt_struct,labStm_struct',[0 0.6])
            else
                plot_identifyabilty_PS(preds_obt_struct,labStm_struct',[0 0.2])
            end
            
        end
        
    end

end