% Make pairwise sequenceness for rest and different phases during task
% GLM approach

function [] = obtain_ClassifierConfusion(include_null,optimise_null,do_normalise,L1_prevBest,...
                                         TS_ELprevBest,do_plot,baseline_correct,Loc_11,which_chan)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==0 

        include_null  = true; % include nulldata when training classifiers
        optimise_null = true; % opmitise amount of nulldata - take as much as needed to minimise correlations between classifiers
        
        % don't change these settings
        do_normalise         = true; % normalise data, should be on        

        temp_smoothing_class = 1;

        L1_prevBest    = 0.006;
        
        TS_ELprevBest  = 20;   
%         TS_RELprevBest = 30;        

        do_plot = true;
%         do_plot = false;        
        
        train_classifier     = false; % do analysis even though file already exists
%         train_classifier     = true; % do analysis even though file already exists        
        
%         baseline_correct = false;
        baseline_correct = true;
        
        Loc_11 = false; % one against all classifiers
%         Loc_11 = true; % one against one classifiers

        which_chan = 'all';
%         which_chan = 'Occ';
%         which_chan = 'Front';
%         which_chan = 'Temp';
%         which_chan = 'Par';
%         which_chan = 'Cen';
        
%         check_accuracy = true; % obtain individual classifier confusion
        check_accuracy = false;

    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP THE MATLAB PATHS and variables
    based = '/Users/Philipp/Dropbox/StimConstrMEG/Results_Scanning/code_MEG_Final';

    scan_result_path    = fullfile(based,'data');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% some relevant variables
    subject = dir(fullfile(scan_result_path,'s*'));
    subject = {subject.name};
    
    idx_trialbeg = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness
    
    fname_betas = 'DataDecode_ALLTS_ClassConf.mat';

    if train_classifier

        %%%%% obtain classifiers %%%%%
        plot_localiser_data = false;

        fprintf('Obtaining element and relation localiser data\n')        
        
        data_class_loc   = nan(4,4,length(subject));

        for idx_sub = 1:length(subject)              

            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, preds_loc,labStm_loc]     = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_localiser.mat',temp_smoothing_class,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);
            
            data_class_loc(:,:,idx_sub)   = plot_identifyabilty_PS(preds_loc,labStm_loc',[0 0.6],false);

            fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))

        end
        
    else

        load(fname_betas)

    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% optional: do stats and plot stuff  
    if do_plot   
         
         plot_mean_identifyabilty_PS(data_class_loc,[0 0.6],'Building block ')
                        
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

end