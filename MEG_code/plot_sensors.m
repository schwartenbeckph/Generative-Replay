% Obtain pairwise sequenceness for inference period
% GLM approach

function [] = plot_sensors(do_localiser,...
                           include_null,optimise_null,...
                           do_normalise,L1_prevBest,TS_ELprevBest,...
                           baseline_correct,Loc_11,which_chan,check_accuracy,plot_localiser_data)
    %% Prelim
    % timing info:
    % Data_inference      [-0.5 3.5]; % epoch end in secs 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==0
        
        which_data = 'inference';
        
        do_normalise = true; % normalise data, should be on
        
        % options important for classifier training
        do_localiser = false; % train classifiers
        include_null  = true; % include nulldata when training classifiers
        optimise_null = true; % opmitise amount of nulldata - take as much as needed to minimise correlations between classifiers    
        L1_prevBest = 0.006;
        TS_ELprevBest  = 20;
%         TS_RELprevBest = 30;   
        baseline_correct = true;        
        Loc_11 = false;
        which_chan = 'all';   
        check_accuracy = false;
        plot_localiser_data = false;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP THE MATLAB PATHS and variables
    % this will be your base directory
    based = '/Users/Philipp/Dropbox/StimConstrMEG/Results_Scanning/code_MEG_Final';

    scan_result_path    = fullfile(based,'data');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% some relevant variables
    subject = dir(fullfile(scan_result_path,'s*'));
    subject = {subject.name};
    
    % these values need to be the same as in the epoching script
    pretrig  = -500;
    fsample  = 100; % after downsampling from 1200Hz
    trloff = round(0.001*pretrig*fsample);             % trloff = how many samples before stim onset
    trlbeg = 0 + trloff;                               % in epoching '0' is replaced by trigger onset, here it's just trloff
    
    idx_trialbeg = abs(trlbeg)+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness 
    
    n_channels = 272;
    
    betas_plot_Sub = zeros(length(subject),n_channels,4);

    for idx_sub = 1:length(subject)   
        
        load(fullfile('data',subject{idx_sub},['Data_',which_data,'.mat'])) % this loads preprocessed 'data' and 'stimlabel'

        %%%%% obtain classifiers %%%%%                
        fprintf('Obtaining element and relation localiser data\n')

        fname_betas = 'Class_data.mat';

        if do_localiser

            [betas_loc, intercepts_loc, corr_betas_loc, amount_null_loc, n_nonZ_chann, betas_loc_11, intercepts_loc_11, idx_chan] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_localiser.mat',1,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);                        
%                     [betas_rel, intercepts_rel, corr_betas_rel, amount_null_rel, n_nonZ_chann_rel, ~, ~, ~, betas_chunk, intercepts_chunk, lab_chunk] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_task.mat',temp_smoothing_class,idx_trialbeg+TS_RELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);
%                     [betas_shape, intercepts_shape, corr_betas_shape, amount_null_shape, n_nonZ_chann_shape, ~, ~, ~, betas_graph, intercepts_graph, lab_graph] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_plan.mat',temp_smoothing_class,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);

        else

            load(fullfile('data',subject{idx_sub},'Class_data.mat'))

        end          

%                 'channel_names','idx_chan_use'

        betas_plot = zeros(length(idx_chan_use),size(betas_loc,2));
        betas_plot(idx_chan_use==1,:) = betas_loc;

        betas_plot_Sub(idx_sub,:,:) = betas_plot;
        
        fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))

    end
    
    zPlotSens_PS(squeeze(mean(mean(betas_plot_Sub,3),1)), 100,[-0.15,0.15])
    zPlotSens_PS(squeeze(mean(betas_plot_Sub(:,:,1),1)), 100,[-0.6,0.6])
    zPlotSens_PS(squeeze(mean(betas_plot_Sub(:,:,2),1)), 100,[-0.6,0.6])
    zPlotSens_PS(squeeze(mean(betas_plot_Sub(:,:,3),1)), 100,[-0.6,0.6])
    zPlotSens_PS(squeeze(mean(betas_plot_Sub(:,:,4),1)), 100,[-0.6,0.6])
  
end