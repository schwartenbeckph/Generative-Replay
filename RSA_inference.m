% Make pairwise sequenceness for rest and different phases during task
% GLM approach

function [] = RSA_inference(which_data,do_RSAana,include_null,optimise_null,do_normalise,L1_prevBest,...
                            TS_ELprevBest,baseline_correct,Loc_11,which_chan)

    % timing info:
    % Data_localiser [-0.5 1.5]; % epoch end in secs  
    % Data_plan      [-0.5 3.5]; % epoch end in secs  
    % Data_task      [-0.5 3.5]; % epoch end in secs 
    % Data_question  [0 1.5]; % epoch end in secs 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==0

        which_data = 'inference';

        include_null     = true; % include nulldata when training classifiers
        optimise_null    = true; % opmitise amount of nulldata - take as much as needed to minimise correlations between classifiers              
        do_normalise     = true; % normalise data, should be on
        L1_prevBest      = 0.006;        
        baseline_correct = true;
        Loc_11           = false;
        which_chan       = 'all';
%         which_chan       = 'Occ';
%         which_chan       = 'Front';
%         which_chan       = 'Temp';
%         which_chan       = 'Par';
%         which_chan       = 'Cen';
        
        do_localiser = false;
%         do_localiser = true;

        do_RSAana     = false;
%         do_RSAana     = true;                 

    end
    
    part_RSA = []; % take all trials
%     part_RSA = 'first_half'; % only first half
%     part_RSA = 'second_half'; % only second half

    n_perm = 5000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP THE MATLAB PATHS and variables
    % this will be your base directory
    based = '/Users/Philipp/Dropbox/StimConstrMEG/Results_Scanning/code_MEG_Final';

    scan_result_path    = fullfile(based,'data');
    behav_result_path   = fullfile(based,'Behav');
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
    
    % find index for onset of stim - used to train classifiers
    idx_trialbeg = abs(trlbeg)+1;  

    trial_max = 6*48;
    
    file_name = 'RSA_inference.mat';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness
    if do_RSAana
        
        Corr_Trials = nan(length(subject),trial_max);
        RTs         = nan(length(subject),trial_max);
        
        r_RSA           = nan(4,4,length(subject));                 
        
        for idx_sub = 1:length(subject)   

            behav_folder   = fullfile(behav_result_path,(subject{idx_sub}));

            %%%%% get biased brick info %%%%%
            [bricks_conn,bricks_rel,bricks_conn_trial,bricks_rel_trial,bricks_conn_Q,bricks_rel_Q,correct_trials_all,rt_all,stims_all,STIMS_all] = mk_bricks_sub(behav_folder); % obtain brick info of subject based on observed stims

            Corr_Trials(idx_sub,1:length(correct_trials_all)) = correct_trials_all;
            RTs(idx_sub,1:length(rt_all))                     = rt_all;
            
            %%%%% Load MEG data %%%%%
            % data: Nchannels x Ntimesteps x Ntrials
            % stimlabel: Ntrials x label         
            load(fullfile(scan_result_path,(subject{idx_sub}),['Data_',which_data,'.mat'])) % this loads 'data' and 'stimlabel'
          
            data = baseline_correct_timeSeries(data,[idx_trialbeg-10,idx_trialbeg-1]); 

            %%%%% obtain classifiers %%%%%
            fprintf('Obtaining element and relation localiser data\n')

            fname_betas = 'Class_data.mat';

            if do_localiser

                [betas_loc, intercepts_loc, corr_betas_loc, amount_null_loc, n_nonZ_chann, betas_loc_11, intercepts_loc_11, idx_chan] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_localiser.mat',1,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);                        
%                 [betas_rel, intercepts_rel, corr_betas_rel, amount_null_rel, n_nonZ_chann_rel, ~, ~, ~, betas_chunk, intercepts_chunk, lab_chunk] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_task.mat',temp_smoothing_class,idx_trialbeg+TS_RELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);
%                 [betas_shape, intercepts_shape, corr_betas_shape, amount_null_shape, n_nonZ_chann_shape, ~, ~, ~, betas_graph, intercepts_graph, lab_graph] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_plan.mat',temp_smoothing_class,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);

            else

                load(fullfile('data',subject{idx_sub},fname_betas))

            end   
            
                fprintf('Running RSA.\n')

                [RSM, EL_RSM, REL_RSM, ~, ~, PIXEL_RSM, SIZE_RSM, ~, ... 
                 ~, ~, ~, ~, ~, ~, ~, ~, ~] = do_RSA_timeseries_inference(fullfile(scan_result_path,(subject{idx_sub})),['Data_',which_data,'.mat'],'all',lab_chunk,STIMS_all,false,part_RSA);              

                EL_RSM    = zscore(EL_RSM);
                REL_RSM   = zscore(REL_RSM);
                PIXEL_RSM = zscore(PIXEL_RSM);
                SIZE_RSM  = zscore(SIZE_RSM);             
                
                smooth_RSA = 15;
                
                for idx_row=1:size(RSM,1)                   
                   RSM(idx_row,:)      = smooth_timeSeries(RSM(idx_row,:),smooth_RSA); 
                end              

                r_RSA(:,:,idx_sub) = corr([EL_RSM REL_RSM PIXEL_RSM SIZE_RSM]);

                if ~exist('RSA_el','var')
                    RSA_el         = nan(length(subject),size(RSM,2),n_perm);
                    RSA_rel        = nan(length(subject),size(RSM,2),n_perm);
                    RSA_pixel      = nan(length(subject),size(RSM,2),n_perm);
                    RSA_size       = nan(length(subject),size(RSM,2),n_perm);
                end                                    
                                                                    
                RSA_el(idx_sub,:,:)    = pick_value_PS(mk_PermutationOLS_PS(RSM,[EL_RSM PIXEL_RSM SIZE_RSM ones(length(EL_RSM),1)],n_perm),...
                                                       {1,1:size(RSM,2),1:n_perm});
                RSA_rel(idx_sub,:,:)   = pick_value_PS(mk_PermutationOLS_PS(RSM,[REL_RSM PIXEL_RSM SIZE_RSM ones(length(EL_RSM),1)],n_perm),...
                                                       {1,1:size(RSM,2),1:n_perm});
                RSA_pixel(idx_sub,:,:) = pick_value_PS(mk_PermutationOLS_PS(RSM,[PIXEL_RSM SIZE_RSM ones(length(EL_RSM),1)],n_perm),...
                                                       {1,1:size(RSM,2),1:n_perm});
                RSA_size(idx_sub,:,:)  = pick_value_PS(mk_PermutationOLS_PS(RSM,[PIXEL_RSM SIZE_RSM ones(length(EL_RSM),1)],n_perm),...
                                                       {2,1:size(RSM,2),1:n_perm});
                                                   
                fprintf('RSA all done.\n')            
            
            fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))

        end
         
         save(file_name,...
             'Corr_Trials',...
             'RTs',...
             'RSA_el','RSA_rel', 'RSA_pixel', 'RSA_size', ...
             'r_RSA', ...            
             '-v7.3')
         
    else
        
        load(file_name)

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plot_timeseries_PS(cat(4,RSA_rel),...
                      'Similarity',[idx_trialbeg-1],{'Building Block in Relation'},...
                      [-500,size(RSA_rel,2)*10-510],'RSA Over Time Planning',true,true,[100 0],[2])


    plot_timeseries_PS(cat(4,RSA_rel,RSA_pixel,RSA_size),...
                      'Similarity',[idx_trialbeg-1],{'Building Block in Relation' 'Shape overlap' 'Size overlap'},...
                      [-500,size(RSA_rel,2)*10-510],'RSA Over Time Planning',true,true,[100 0],[2 3 5 6])

end