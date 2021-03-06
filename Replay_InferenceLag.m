% Make pairwise sequenceness for rest and different phases during task
% GLM approach

function [] = Replay_InferenceLag(which_data,include_null,optimise_null,do_normalise,L1_prevBest,...
                                  TS_ELprevBest,do_noise,control_oscillation,tw,baseline_correct,Loc_11,which_chan)

    % timing info:
    % Data_localiser [-0.5 1.5]; % epoch end in secs  
    % Data_plan      [-0.5 3.5]; % epoch end in secs  
    % Data_task      [-0.5 3.5]; % epoch end in secs 
    % Data_question  [0 1.5]; % epoch end in secs 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==0
 
        which_data = 'inference';
        
        do_normalise = true; % normalise data, should be on
        
        do_noise = false; % turn on if do analysis with simulated noise
        
        % we can control for oscillations, but not here since not resting
        % data
%         control_oscillation = 'ConAlpha';
%         control_oscillation = 'ConTheta';
        control_oscillation = 'ConNo';         

        do_replayAna = false; % do analysis even though file already exists
%         do_replayAna = true; % do analysis even though file already exists
        
        % options important for classifier training
%         do_localiser = true; % train classifiers
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
    
    combine_sess = 6; % how many sessions to combine when looking at task data, 1 == don't combine

    n_perm = 3;
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

    nstates = 4; % 4 elements
    
    % these values need to be the same as in the epoching script
    pretrig  = -500;
    fsample  = 100; % after downsampling from 1200Hz
    trloff = round(0.001*pretrig*fsample);             % trloff = how many samples before stim onset
    trlbeg = 0 + trloff;                               % in epoching '0' is replaced by trigger onset, here it's just trloff
    
    % find index for onset of stim - used to train classifiers
    idx_trialbeg = abs(trlbeg)+1;    
    
    maxLag = 50; % evaluate sequenceness up to 500ms
          
    trial_max = 6*48;
    
    % Filename of pre-analysed sequence data
    file_name = 'Replay_InferenceLag.mat';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness    
%     n_shuffles = 1;
    
    if do_replayAna
        
        s_ElEl       = nan(length(subject),maxLag,trial_max,n_perm,3);
        s_ElEl_wrong = nan(length(subject),maxLag,trial_max,n_perm,4);        
        
        Corr_Trials = nan(length(subject),trial_max);
        RTs         = nan(length(subject),trial_max);        
        
        betasnbins64_allSub = nan(maxLag,nstates^2,length(subject),trial_max);
        preds_all_allSub = nan(length(subject),length(tw),nstates);

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
            
            if do_noise
                data = rand(size(data));
            end            
            
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
            
            % reduce data to some sensors
            data = data(idx_chan,:,:);           

            %%%%% loop over trial %%%%%
            use_trials = 1:length(rt_all);
            use_trials = use_trials(~isnan(squeeze(data(1,idx_trialbeg,:))));
            
            if ~exist('ERF_subj','var')
                ERF_subj = nan(length(subject),size(data,2));
            end
            
            ERF_subj(idx_sub,:) = squeeze(nanmean(nanmean(data(:,1:size(ERF_subj,2),:),3),1));
            if do_normalise
                ERF_subj(idx_sub,:) = mk_normalise(ERF_subj(idx_sub,:));
            end

            for idx_trial = use_trials                                       
                       
               if isempty(tw)
                   X = data(:,idx_trialbeg:end,idx_trial);
               elseif length(tw)==1
                   X = data(:,idx_trialbeg+tw:end,idx_trial);
               else
                   X = data(:,idx_trialbeg+tw,idx_trial);
               end

               X = nanmean(X,3)'; % concatenate instead of mean?
                   
               if ~exist('preds_present')
                   preds_present = nan(size(X,1),length(subject),trial_max);
                   preds_absent = nan(size(X,1),length(subject),trial_max);
               end

               X(isnan(X(:,1)),:) = []; 
               
               if ~isempty(X)

                   X = mk_normalise(X);

                   %%%%% make prediction for state at any given time point %%%%%
                   % this is time-point x class                       
                   preds_el = 1./(1+exp(-(X*betas_loc+ repmat(intercepts_loc, [size(X,1) 1]))));

                   preds_all = preds_el;

                   preds_all_allSub(idx_sub,1:size(preds_all,1),:) = preds_all;
                   
                   el_present = unique(bricks_rel_trial(idx_trial,:),'stable');
                   el_absent  = setdiff([1:4],el_present);
                   
                   preds_present(:,idx_sub,idx_trial) = mean(preds_el(:,el_present),2);
                   preds_absent(:,idx_sub,idx_trial)  = mean(preds_el(:,el_absent),2);                    
                   
                   nbins = maxLag+1;

                   %%%%% take the predictions for each of the four classes and shift them by iLag=1:maxLag %%%%% 
                   % this is time-point x (#classes * maxLag), columns 1:60 are time shifted class 1, then class 2, etc...
                   warning off

                   dm = [];
                   for kk=1:nstates % do for all classes 
                       temp = toeplitz(preds_all(:,kk),zeros(nbins,1)); % this is time-point x maxLag+1 (includes zero-lag predicitons)
                       temp = temp(:,2:end); % this is time-point x maxLag (exclude zero-lag predicitons)
                       dm   = [dm temp]; 
                   end
                                      
                   warning on
                   % now we have a big matrix that contains predicted classes at all possible time lags

                   betas = nan(nstates*maxLag, nstates);

                   %%%%% "First-level GLM": control for other lags (alpha) %%%%% 
                   if strcmp(control_oscillation,'ConAlpha')
                       bins = 10; % control for alpha 100ms cycles
                   elseif strcmp(control_oscillation,'ConTheta')
                       bins = 25; % theta
                   elseif strcmp(control_oscillation,'ConNo')
                       bins = maxLag; % no control
                   end

                   % does the time-lagged data 'predict itself' - i.e. 
                   for ilag=1:bins
                       temp_zinds          = (1:bins:nstates*maxLag) + ilag - 1; % now just cut out different columns of dm <=> different time lags
                       temp                = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*preds_all; % obtain predictors for different time lags for observed pattern
                       betas(temp_zinds,:) = temp(1:end-1,:);    
                   end
                   % 'betas' tells us: how good can we predict the four
                   % classes with time-lagged data of the four classes
                   % where:
                   % rows 1:maxLag is how well do we predict the four
                   % classes with time lagged data from class 1
                   % rows maxLag+1:maxLag*2 is how well do we predict the four
                   % classes with time lagged data from class 2
                   % ... same for class 3 and class 4

                   %%%%% now you have a beta for every lag*class for the four classes, i.e. a transition pattern over time! %%%%%
                   betasnbins64 = reshape(betas,[maxLag nstates^2]); % reshape to turn this into transition pattern per time lag
                   
                   betasnbins64_allSub(:,:,idx_sub,idx_trial) = betasnbins64;
            
                   % control transitions:
                   % mean
                   T_mean = ones(nstates); 
                   
                   % self transitions
                   T_self = eye(nstates);
                   
                   T_self_El          = zeros(nstates); 
                   T_self_El(1:4,1:4) = T_self(1:4,1:4);
                   
                   T_self_all = T_self_El(:);
                   
                   % decide which mean to use
                   T_mean_use = T_mean(:);
                       
                   [T_ElEl,~,~,T_EE_wrong,~,~,~,~,~] = mk_experienced_transitions(bricks_conn_trial(idx_trial,:),bricks_rel_trial(idx_trial,:),nstates);                       

                   % CHECK whether element 4 is always the one that is always there
                   BB_alwaysPrsnt = find([all(sum(ismember(bricks_conn_trial,1),2)) all(sum(ismember(bricks_conn_trial,3),2)) ...
                                          all(sum(ismember(bricks_conn_trial,3),2)) all(sum(ismember(bricks_conn_trial,4),2))]); 

                   BB_vary = setdiff(1:4,BB_alwaysPrsnt);

                   if idx_trial==1
                        fprintf('Always present element is %d.\n',BB_alwaysPrsnt)
                   end

                   T_ElEl_stableTovary    = T_ElEl; T_ElEl_stableTovary(BB_vary,:) = 0;
                   T_ElEl_varyTostable    = T_ElEl; T_ElEl_varyTostable(:,BB_vary) = 0;
                   % this may sometimes be all zeros:
                   T_ElEl_varyTovary      = T_ElEl; T_ElEl_varyTovary(BB_alwaysPrsnt,:)     = 0; T_ElEl_varyTovary(:,BB_alwaysPrsnt) = 0;                                                                      
                                       
                   T_ElEl_stableTovary_alt = zeros(size(T_ElEl_stableTovary,1),size(T_ElEl_stableTovary,2),2);
                   idx_stableTovary        = find(T_ElEl_stableTovary);
                   idx_stableTovary_alt    = [4 8 12];
                   idx_stableTovary_alt    = setdiff(idx_stableTovary_alt,idx_stableTovary);
                   if length(idx_stableTovary)==2
                       T_ElEl_stableTovary_alt([idx_stableTovary(1) idx_stableTovary_alt]) = 1;
                       T_ElEl_stableTovary_alt([idx_stableTovary(2) idx_stableTovary_alt]+length(T_ElEl_stableTovary(:))) = 1;
                   else
                       T_ElEl_stableTovary_alt([idx_stableTovary_alt(1)]) = 1;
                       T_ElEl_stableTovary_alt([idx_stableTovary_alt(2)]+length(T_ElEl_stableTovary(:))) = 1;
                   end

%                    figure,imagesc(T_ElEl_stableTovary)
%                    figure,imagesc(T_ElEl_stableTovary_alt(:,:,1))
%                    figure,imagesc(T_ElEl_stableTovary_alt(:,:,2))

                   T_ElEl_varyTostable_alt = zeros(size(T_ElEl_varyTostable,1),size(T_ElEl_varyTostable,2),2); 
                   idx_varyTostable        = find(T_ElEl_varyTostable);
                   idx_varyTostable_alt    = [13 14 15];
                   idx_varyTostable_alt    = setdiff(idx_varyTostable_alt,idx_varyTostable);
                   if length(idx_varyTostable)==2
                       T_ElEl_varyTostable_alt([idx_varyTostable(1) idx_varyTostable_alt]) = 1;
                       T_ElEl_varyTostable_alt([idx_varyTostable(2) idx_varyTostable_alt]+length(T_ElEl_varyTostable(:))) = 1;
                   else
                       T_ElEl_varyTostable_alt([idx_varyTostable_alt(1)]) = 1;
                       T_ElEl_varyTostable_alt([idx_varyTostable_alt(2)]+length(T_ElEl_varyTostable(:))) = 1;
                   end       

%                    figure,imagesc(T_ElEl_varyTostable)
%                    figure,imagesc(T_ElEl_varyTostable_alt(:,:,1))
%                    figure,imagesc(T_ElEl_varyTostable_alt(:,:,2))


                   if any(T_ElEl_varyTovary(:))==1

                       T_ElEl_varyTovary_alt = zeros(size(T_ElEl_varyTovary,1),size(T_ElEl_varyTovary,2),2); 
                       idx_varyTovary        = find(T_ElEl_varyTovary);
                       idx_varyTovary_alt    = [2 5;
                                                3 9;
                                                7 10];
                       idx_varyTovary_alt    = idx_varyTovary_alt(~ismember(idx_varyTovary_alt,idx_varyTovary','rows'),:);
                       T_ElEl_varyTovary_alt(idx_varyTovary_alt(1,:)) = 1;
                       T_ElEl_varyTovary_alt(idx_varyTovary_alt(2,:)+length(T_ElEl_varyTovary(:))) = 1;

%                        figure,imagesc(T_ElEl_varyTovary)
%                        figure,imagesc(T_ElEl_varyTovary_alt(:,:,1))
%                        figure,imagesc(T_ElEl_varyTovary_alt(:,:,2))

                   else

                       T_ElEl_varyTovary_alt = [];

                   end
                   
                   T_EE_wrong_stableToabsent    = T_EE_wrong; T_EE_wrong_stableToabsent(BB_vary,:) = 0;
                   % this may sometimes be all zeros:
                   T_EE_wrong_varyTovary      = T_EE_wrong; T_EE_wrong_varyTovary(BB_alwaysPrsnt,:)     = 0; T_EE_wrong_varyTovary(:,BB_alwaysPrsnt) = 0;
                   T_EE_wrong_absentTovary    = T_EE_wrong_varyTovary; T_EE_wrong_absentTovary(el_present,:) = 0;
                   T_EE_wrong_varyToabsent    = T_EE_wrong_varyTovary; T_EE_wrong_varyToabsent(el_absent,:) = 0;

                   T_EE_wrong_absentTostable    = T_EE_wrong; T_EE_wrong_absentTostable(:,BB_vary) = 0;

                                   
                   if all(T_ElEl_varyTovary(:)==0)
                       

                       bbb = pinv([T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl(idx_sub,:,idx_trial,1,2)    = bbb(1,:);
                       s_ElEl(idx_sub,:,idx_trial,1,3)    = bbb(2,:); 

                       bbb = pinv([squash(T_ElEl_varyTostable_alt(:,:,1)) T_ElEl_stableTovary(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,2,2)    = bbb(1,:);
                       bbb = pinv([squash(T_ElEl_varyTostable_alt(:,:,2)) T_ElEl_stableTovary(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,3,2)    = bbb(1,:);

                       bbb = pinv([T_ElEl_varyTostable(:) squash(T_ElEl_stableTovary_alt(:,:,1)) ...
                                   T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,2,3)    = bbb(2,:);
                       bbb = pinv([T_ElEl_varyTostable(:) squash(T_ElEl_stableTovary_alt(:,:,2)) ...
                                   T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,3,3)    = bbb(2,:);

                       bbb = pinv([T_EE_wrong_absentTovary(:) T_EE_wrong_varyToabsent(:) ...
                                   T_EE_wrong_absentTostable(:) T_EE_wrong_stableToabsent(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl_wrong(idx_sub,:,idx_trial,1,1)    = bbb(1,:);
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,2)    = bbb(2,:); 
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,3)    = bbb(3,:);
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,4)    = bbb(4,:);

                   else                                           

                       bbb = pinv([T_ElEl_varyTovary(:) T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl(idx_sub,:,idx_trial,1,1)    = bbb(1,:);
                       s_ElEl(idx_sub,:,idx_trial,1,2)    = bbb(2,:);
                       s_ElEl(idx_sub,:,idx_trial,1,3)    = bbb(3,:);

                       bbb = pinv([squash(T_ElEl_varyTovary_alt(:,:,1)) T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,2,1)    = bbb(1,:);
                       bbb = pinv([squash(T_ElEl_varyTovary_alt(:,:,2)) T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,3,1)    = bbb(1,:);

                       bbb = pinv([T_ElEl_varyTovary(:) squash(T_ElEl_varyTostable_alt(:,:,1)) T_ElEl_stableTovary(:) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,2,2)    = bbb(2,:);
                       bbb = pinv([T_ElEl_varyTovary(:) squash(T_ElEl_varyTostable_alt(:,:,2)) T_ElEl_stableTovary(:) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,3,2)    = bbb(2,:);

                       bbb = pinv([T_ElEl_varyTovary(:) T_ElEl_varyTostable(:) squash(T_ElEl_stableTovary_alt(:,:,1)) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,2,3)    = bbb(3,:);
                       bbb = pinv([T_ElEl_varyTovary(:) T_ElEl_varyTostable(:) squash(T_ElEl_stableTovary_alt(:,:,2)) ...
                               T_mean_use T_self_all])*(betasnbins64');
                       s_ElEl(idx_sub,:,idx_trial,3,3)    = bbb(3,:);

                       bbb = pinv([T_EE_wrong_absentTovary(:) T_EE_wrong_varyToabsent(:) ...
                                   T_EE_wrong_absentTostable(:) T_EE_wrong_stableToabsent(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl_wrong(idx_sub,:,idx_trial,1,1)    = bbb(1,:);
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,2)    = bbb(2,:); 
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,3)    = bbb(3,:);
                       s_ElEl_wrong(idx_sub,:,idx_trial,1,4)    = bbb(4,:);

                   end
               
               end               

            end
            
            fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))
            
        end
         
    else
        
        load(file_name)

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% optional: do stats and plot stuff    
    % plotting colours
    col   = {[0, 0.4470, 0.7410], ...       % blue
             [0.4940, 0.1840, 0.5560], ...  % purple
             [0.4660, 0.6740, 0.1880], ...  % green
             [0.9350, 0.1780, 0.2840], ...  % red
             [0.3010, 0.7450, 0.9330], ...  % cyan
             [0.9290, 0.6940, 0.1250], ...  % yellow
             [0, 0, 0], ...                 % black
             [0.6, 0.6, 0.6]};              % grey 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get means of sessions during task
    s_ElEl       = mk_MeanSess(s_ElEl); 

    s_ElEl       = mk_CombineSess(s_ElEl,combine_sess);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now plot stuff   
    plot_what = 'mean';       
    
    plot_diff_sequenceness(cat(5,nanmean(s_ElEl(:,:,:,:,1:2),5),s_ElEl(:,:,:,:,3)),[],...
                           ['Building Block to Building Block Seqs ',which_data],10000,true,false,true,...
                           {'Not always present -> ' ...
                            'Always present ->'},plot_what,true,false,[],[],false,false,'plan',false,[7 4])                                                        

end