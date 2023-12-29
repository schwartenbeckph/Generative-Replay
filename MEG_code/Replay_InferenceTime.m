% Obtain pairwise sequenceness for inference period
% GLM approach

function [] = Replay_InferenceTime(which_data,do_replayAna,do_localiser,...
                                   include_null,optimise_null,...
                                   do_normalise,L1_prevBest,TS_ELprevBest,...
                                   do_noise,baseline_correct,Loc_11,which_chan,check_accuracy,plot_localiser_data)
    %% Prelim
    % timing info:
    % Data_inference      [-0.5 3.5]; % epoch end in secs 
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
        do_localiser = false; % train classifiers
        include_null  = true; % include nulldata when training classifiers
        optimise_null = true; % opmitise amount of nulldata - take as much as needed to minimise correlations between classifiers    
        L1_prevBest = 0.006;
        TS_ELprevBest  = 20;  
        baseline_correct = true;        
        Loc_11 = false;
        which_chan = 'all';   
        check_accuracy = false;
        plot_localiser_data = false;

    end    

    % define time windows
    start_tw = -50;
    end_tw   = 350;

    idx_incr = 1; % this reproduces main finding in figure 6
%     idx_incr = 50; % this reproduces supplementary figure 5
%     idx_incr = 25; % this reproduces supplementary figure 5

    window_size = 100;

    idx_tw = 1;
    idx_tw_end = start_tw+window_size;
    while idx_tw_end<=end_tw
        tw_all{idx_tw} = start_tw:idx_tw_end;

        idx_tw = idx_tw+1;

        start_tw = start_tw+idx_incr;

        idx_tw_end = start_tw+window_size;
    end

    combine_sess = 6; % how many sessions to combine when looking at task data, 1 == don't combine    
    
    n_perm = 1; % no permutations of sequences here, but see Replay_TimeLags.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SETUP THE MATLAB PATHS and variables
    % this will be your base directory
    based = '';

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
    
    idx_trialbeg = abs(trlbeg)+1;
    
    maxLag = 20; % evaluate sequenceness up to 200ms 
    
    trial_max = 6*48;
    
    % Filename of pre-analysed sequence data
    if idx_incr == 1
        file_name = 'Replay_InferenceTime.mat'; % sliding window
    elseif idx_incr == 50
        file_name = 'Replay_InferenceTime_SepIntervals.mat'; % separate intervals
    elseif idx_incr == 25
        file_name = 'idx_incr25_Replay_InferenceTime_SepIntervals.mat'; % separate intervals
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness    
    distances_seq = nan(length(subject),length(tw_all));

    vary_BB_all     = nan(length(subject),maxLag,length(tw_all));
    Distvary_BB_all = nan(length(subject),maxLag,length(tw_all));
    stable_BB_all   = nan(length(subject),maxLag,length(tw_all));    
    BB_tw_all       = nan(length(subject),maxLag,length(tw_all),4);

    vary_BB_wrong_all   = nan(length(subject),maxLag,length(tw_all)); 
    stable_BB_wrong_all = nan(length(subject),maxLag,length(tw_all)); 
    BB_wrong_tw_all     = nan(length(subject),maxLag,length(tw_all),4); 
    
    if do_replayAna
    
        for idx_timeW = 1:length(tw_all)

            tw = tw_all{idx_timeW};      

            s_ElEl       = nan(length(subject),maxLag,trial_max,n_perm,4);
            s_ElEl_wrong = nan(length(subject),maxLag,trial_max,n_perm,4);

            Corr_Trials = nan(length(subject),trial_max);
            RTs         = nan(length(subject),trial_max);        
   
            betasnbins64_allSub = nan(maxLag,nstates^2,length(subject),trial_max);

            preds_all_allSub = nan(length(subject),length(tw),nstates);

            for idx_sub = 1:length(subject)   

                behav_folder   = fullfile(behav_result_path,(subject{idx_sub}));

                %%%%% get biased brick info %%%%%
                [~,~,bricks_conn_trial,bricks_rel_trial,bricks_conn_Q,bricks_rel_Q,correct_trials_all,rt_all,~,~] = mk_bricks_sub(behav_folder); % obtain brick info of subject based on observed stims

                Corr_Trials(idx_sub,1:length(correct_trials_all)) = correct_trials_all;
                RTs(idx_sub,1:length(rt_all))                     = rt_all;            

                %%%%% Load MEG data %%%%%
                % data: Nchannels x Ntimesteps x Ntrials
                % stimlabel: Ntrials x label         
                load(fullfile('data',subject{idx_sub},['Data_',which_data,'.mat'])) % this loads preprocessed 'data' and 'stimlabel'
                
                if do_noise
                    data = rand(size(data));
                end            

                data = baseline_correct_timeSeries(data,[idx_trialbeg-10,idx_trialbeg-1]); 

                %%%%% obtain classifiers %%%%%                
                fprintf('Obtaining element and relation localiser data\n')

                fname_betas = 'Class_data.mat';
                
                if do_localiser

                    [betas_loc, intercepts_loc, corr_betas_loc, amount_null_loc, n_nonZ_chann, betas_loc_11, intercepts_loc_11, idx_chan] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_localiser.mat',1,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);                        
%                     [betas_rel, intercepts_rel, corr_betas_rel, amount_null_rel, n_nonZ_chann_rel, ~, ~, ~, betas_chunk, intercepts_chunk, lab_chunk] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_task.mat',temp_smoothing_class,idx_trialbeg+TS_RELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);
%                     [betas_shape, intercepts_shape, corr_betas_shape, amount_null_shape, n_nonZ_chann_shape, ~, ~, ~, betas_graph, intercepts_graph, lab_graph] = obtain_betas(fullfile(scan_result_path,(subject{idx_sub})),'Data_plan.mat',temp_smoothing_class,idx_trialbeg+TS_ELprevBest,include_null,optimise_null,L1_prevBest,do_normalise,[],[],plot_localiser_data,baseline_correct,Loc_11,which_chan,check_accuracy);

                else
                    
                    load(fullfile('data',subject{idx_sub},fname_betas))

                end                          

                % reduce data to some sensors
                data = data(idx_chan,:,:);           

                %%%%% loop over trial %%%%%
                use_trials = 1:length(rt_all);
                use_trials = use_trials(~isnan(squeeze(data(1,idx_trialbeg,:))));

                % plot ERF:
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

                   X = nanmean(X,3)';            

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

                   end

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
                   T_self       = eye(nstates);

                   T_self_El          = zeros(nstates); 
                   T_self_El(1:4,1:4) = T_self(1:4,1:4);

                   T_self_all = [T_self_El(:)];

                   T_mean_use = T_mean(:);

                   % obtain transition probs between present building blocks and
                   % absent ('wrong') building blocks
                   [T_ElEl,~,~,T_EE_wrong,~,~,~,~,~] = mk_experienced_transitions(bricks_conn_trial(idx_trial,:),bricks_rel_trial(idx_trial,:),nstates);

                   BB_alwaysPrsnt = find([all(sum(ismember(bricks_conn_trial,1),2)) all(sum(ismember(bricks_conn_trial,3),2)) ...
                                          all(sum(ismember(bricks_conn_trial,3),2)) all(sum(ismember(bricks_conn_trial,4),2))]); 

                   BB_vary = setdiff(1:4,BB_alwaysPrsnt);

                   if idx_trial==1
                        fprintf('Always present element is %d.\n',BB_alwaysPrsnt) % building block 4 is always present
                   end

                   % now split those building blocks up
                   T_ElEl_stableTovary    = T_ElEl; T_ElEl_stableTovary(BB_vary,:) = 0; % stable -> present
                   T_ElEl_varyTostable    = T_ElEl; T_ElEl_varyTostable(:,BB_vary) = 0; % present -> stable
                   % this may sometimes be all zeros:
                   T_ElEl_varyTovary      = T_ElEl; 
                   T_ElEl_varyTovary(BB_alwaysPrsnt,:) = 0; 
                   T_ElEl_varyTovary(:,BB_alwaysPrsnt) = 0; % present <-> present  
                   
                   T_ElEl_distantvaryTostable = zeros(size(T_ElEl));
                   
                   if any(T_ElEl_varyTovary(:))==1
                       distantTostable = setdiff(find(sum(T_ElEl_varyTovary,2)),find(sum(T_ElEl_varyTostable,2)));                       
                       T_ElEl_distantvaryTostable(distantTostable,BB_alwaysPrsnt) = 1; % distant (unconnected) present -> stable
                   end

                   % same for absent building block
                   T_EE_wrong_stableToabsent    = T_EE_wrong; T_EE_wrong_stableToabsent(BB_vary,:) = 0; % stable -> absent
                   
                   % this may sometimes be all zeros:
                   T_EE_wrong_varyTovary      = T_EE_wrong; T_EE_wrong_varyTovary(BB_alwaysPrsnt,:)     = 0; T_EE_wrong_varyTovary(:,BB_alwaysPrsnt) = 0; % absent <-> present
                   T_EE_wrong_absentTovary    = T_EE_wrong_varyTovary; T_EE_wrong_absentTovary(el_present,:) = 0; % absent -> present
                   T_EE_wrong_varyToabsent    = T_EE_wrong_varyTovary; T_EE_wrong_varyToabsent(el_absent,:) = 0; % present -> absent

                   T_EE_wrong_absentTostable    = T_EE_wrong; T_EE_wrong_absentTostable(:,BB_vary) = 0; % absent -> stable
                   
                   
                   % if stable building block is in the middle
                   if all(T_ElEl_varyTovary(:)==0)

                       bbb = pinv([T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl(idx_sub,:,idx_trial,1,2)    = bbb(1,:);
                       s_ElEl(idx_sub,:,idx_trial,1,3)    = bbb(2,:);


                   else                                           

                       bbb = pinv([T_ElEl_varyTovary(:) T_ElEl_varyTostable(:) T_ElEl_stableTovary(:) T_ElEl_distantvaryTostable(:) ...
                                   T_mean_use T_self_all])*(betasnbins64');

                       s_ElEl(idx_sub,:,idx_trial,1,1)    = bbb(1,:);
                       s_ElEl(idx_sub,:,idx_trial,1,2)    = bbb(2,:);
                       s_ElEl(idx_sub,:,idx_trial,1,3)    = bbb(3,:);
                       s_ElEl(idx_sub,:,idx_trial,1,4)    = bbb(4,:);

                   end
                   
                   bbb = pinv([T_EE_wrong_absentTovary(:) T_EE_wrong_varyToabsent(:) ...
                   T_EE_wrong_absentTostable(:) T_EE_wrong_stableToabsent(:) ...
                   T_mean_use T_self_all])*(betasnbins64');

                   s_ElEl_wrong(idx_sub,:,idx_trial,1,1)    = bbb(1,:); % absent -> vary
                   s_ElEl_wrong(idx_sub,:,idx_trial,1,2)    = bbb(2,:); % vary -> absent
                   s_ElEl_wrong(idx_sub,:,idx_trial,1,3)    = bbb(3,:); % absent -> stable
                   s_ElEl_wrong(idx_sub,:,idx_trial,1,4)    = bbb(4,:); % stable -> absent

                end               

            end

            fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))
            
            s_ElEl       = mk_MeanSess(s_ElEl);            
            s_ElEl_wrong = mk_MeanSess(s_ElEl_wrong);
            
            s_ElEl       = mk_CombineSess(s_ElEl,combine_sess);
            s_ElEl_wrong = mk_CombineSess(s_ElEl_wrong,combine_sess);
            
            vary_BB     = squeeze(nanmean(s_ElEl(:,:,:,:,1:2),5));
            stable_BB   = squeeze(s_ElEl(:,:,:,:,3));
            Distvary_BB = squeeze(s_ElEl(:,:,:,:,4));
            
            vary_BB_wrong   = squeeze(nanmean(s_ElEl_wrong(:,:,:,:,1:2),5));
            stable_BB_wrong = squeeze(nanmean(s_ElEl_wrong(:,:,:,:,3:4),5));
            
            distances_seq(:,idx_timeW,:) = nanmean(vary_BB(:,1:20,:),2)-nanmean(stable_BB(:,1:20,:),2);            
            
            vary_BB_all(:,:,idx_timeW,:)     = vary_BB; 
            stable_BB_all(:,:,idx_timeW,:)   = stable_BB; 
            Distvary_BB_all(:,:,idx_timeW,:) = Distvary_BB; 
            
            BB_tw_all(:,:,idx_timeW,:,:) = squeeze(s_ElEl);
            
            vary_BB_wrong_all(:,:,idx_timeW,:)   = vary_BB_wrong; 
            stable_BB_wrong_all(:,:,idx_timeW,:) = stable_BB_wrong; 
            
            BB_wrong_tw_all(:,:,idx_timeW,:,:) = squeeze(s_ElEl_wrong);   
            
            fprintf('########Done with tw %d of %d.########\n',idx_timeW,length(tw_all))
        
        end
    
    else

        load(file_name)

    end         

    Abs_var    = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,1),4),2)); % absent -> vary
    Var_abs    = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,2),4),2)); % vary -> absent
    Abs_stable = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,3),4),2)); % absent -> stable 
    
    Var_var    = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,1),4),2)); % vary -> vary
    Var_stable = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,2),4),2)); % vary -> stable
    Stable_var = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,3),4),2)); % stable -> vary
    DistVar_stable = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,4),4),2)); % stable -> vary        

    if size(distances_seq,2) == 7
        
        [h_dist, ~] = ttest(distances_seq);
    
        [h_VarVar, ~]  = ttest(Var_var);
        [h_VarStab, ~] = ttest(Var_stable);
        [h_StabVar, ~] = ttest(Stable_var);
        [h_DistVarStab, ~] = ttest(DistVar_stable);
        
        Present_data = cat(3,Var_var,Var_stable,Stable_var,DistVar_stable);
    
        mk_PlotSotErrbar_Group_PS(distances_seq,true,'To-Be-Inferred - Always Present','Time-Interval','Difference To-Be-Inferred - Always Present Different Time-Intervals',...
                                  {'-500ms - 500ms' '0ms - 1000ms' '500ms - 1500ms' '1000ms - 2000ms' '1500ms - 2500ms' '2000ms - 2500ms' '2500ms - 3500ms'},...
                                  h_dist,1,[],[-0.05 0.05],[0 8],[])
                              
        mk_PlotSotErrbar_Group_PS(Present_data,true,' ',...
                                  'Time-Interval','Difference To-Be-Inferred - Always Present Different Time-Intervals',...
                                  {'-500ms - 500ms' '0ms - 1000ms' '500ms - 1500ms' '1000ms - 2000ms' '1500ms - 2500ms' '2000ms - 2500ms' '2500ms - 3500ms'},...
                                  cat(3,h_VarVar,h_VarStab,h_StabVar,h_DistVarStab),[2 3 4 6],{'vary -> vary' 'vary -> stable' 'stable -> vary' 'dist vary -> stable'},...
                                  [-0.05 0.05],[0 8],[])      
    elseif size(distances_seq,2) == 13
        
        [h_VarVar, ~]  = ttest(Var_var);
        [h_VarStab, ~] = ttest(Var_stable);
        [h_StabVar, ~] = ttest(Stable_var);
        [h_DistVarStab, ~] = ttest(DistVar_stable);
        
        Present_data = cat(3,Var_var,Var_stable,Stable_var,DistVar_stable);
                              
        mk_PlotSotErrbar_Group_PS(Present_data,true,' ',...
                                  'Time-Interval','Difference To-Be-Inferred - Always Present Different Time-Intervals',...
                                  {'-500ms - 500ms' '-250ms - 750ms' '0ms - 1000ms' '250ms - 1250ms' '500ms - 1500ms' '750ms - 1750ms' '1000ms - 2000ms' '1250ms - 2250ms' '1500ms - 2500ms' '1750ms - 2750ms'  '2000ms - 3000ms'  '2250ms - 3250ms'  '2500ms - 3500ms'},...
                                  cat(3,h_VarVar,h_VarStab,h_StabVar,h_DistVarStab),[2 3 4 6],{'vary -> vary' 'vary -> stable' 'stable -> vary' 'dist vary -> stable'},...
                                  [-0.05 0.05],[0 14],[]) 
                              
    else
        
        % plot difference
%         [~, ~, ~, sig_t_maxClus] = sign_flip_AcrossTime_PS(distances_seq,10000,[],[100,0],true,0.05,false);
    %     Note: the above is noisy based on 10000 shuffles, run below to get more reliable estimate
        [~, ~, ~, sig_t_maxClus] = sign_flip_AcrossTime_PS(distances_seq,200000,[],[100,0],true,0.05,false); 

        [sig_where] = mk_clusterStats(distances_seq,sig_t_maxClus,true,false,size(distances_seq,1)-1,0.05,true,false,9,'Seqs vary vs. stable BBs different intervals',[],[-0.05,0.04]);

        time = -50:350;
        disp(sig_where)
        time(sig_where{1})
        
        % plot individual variables                          
        [~, ~, ~, sig_t_maxClus_VarVar]      = sign_flip_AcrossTime_PS(Var_var,10000,[],[100,0],true,0.05,false);            
        [~, ~, ~, sig_t_maxClus_VarStab]     = sign_flip_AcrossTime_PS(Var_stable,10000,[],[100,0],true,0.05,false);   
        [~, ~, ~, sig_t_maxClus_StableVar]   = sign_flip_AcrossTime_PS(Stable_var,10000,[],[100,0],true,0.05,false); 
        [~, ~, ~, sig_t_maxClus_AbsStable]   = sign_flip_AcrossTime_PS(Abs_stable,10000,[],[100,0],true,0.05,false);
        [~, ~, ~, sig_t_maxClus_DistVarStab] = sign_flip_AcrossTime_PS(DistVar_stable,10000,[],[100,0],true,0.05,false);

        [sig_where] =  mk_clusterStats(cat(3,Var_var,Stable_var,Abs_stable,DistVar_stable,Var_stable),...
                                       cat(3,sig_t_maxClus_VarVar,sig_t_maxClus_StableVar,sig_t_maxClus_AbsStable,sig_t_maxClus_DistVarStab,sig_t_maxClus_VarStab),...
                                       true,false,size(Stable_var,1)-1,0.05,true,false,[2 4 5 6 3],'Individual Sequences',...
                                       {'Present - present','Stable - present','Abssent - stable','Distant Present - stable','Presend - stable'},...
                                       [-0.07,0.07]);  
                                   
        disp(sig_where)
        time(sig_where{1}{3})
        time(sig_where{1}{4})
        time(sig_where{1}{5})

        % plot absent stuff           
        AbsVarVarAbs = mean(cat(3,Abs_var,Var_abs),3);
        [~, ~, ~, sig_t_maxClus_AbsVarVarAbs]    = sign_flip_AcrossTime_PS(AbsVarVarAbs,10000,[],[100,0],true,0.05,false);  

        difference_VarVar_AbsVarVarAbs = Var_var - mean(cat(3,Abs_var,Var_abs),3);
        [~, ~, ~, sig_t_maxClus_VarVar_AbsVarVarAbs]    = sign_flip_AcrossTime_PS(difference_VarVar_AbsVarVarAbs,10000,[],[100,0],true,0.05,false);   


        [sig_where] =  mk_clusterStats(cat(3,difference_VarVar_AbsVarVarAbs,Var_var,AbsVarVarAbs),...
                                       cat(3,sig_t_maxClus_VarVar_AbsVarVarAbs,sig_t_maxClus_VarVar,sig_t_maxClus_AbsVarVarAbs),...
                                       true,false,size(Stable_var,1)-1,0.05,true,false,[7 2 8],'Individual Seqs',...
                                       {'Difference','Present - present','Present - absent'},...
                                       [-0.03,0.04]);  
                                   
        disp(sig_where)
        time(sig_where{1}{1})

        % plot all sig stuff
        [sig_where] =  mk_clusterStats(cat(3,Abs_stable,DistVar_stable,Var_stable,difference_VarVar_AbsVarVarAbs),...
                                       cat(3,sig_t_maxClus_AbsStable,sig_t_maxClus_DistVarStab,sig_t_maxClus_VarStab,sig_t_maxClus_VarVar_AbsVarVarAbs),...
                                       true,false,size(Stable_var,1)-1,0.05,true,false,[5 6 3,7],'Individual Seqs',...
                                       {'Absent - stable','Distant Present - stable','Presend - stable','Difference'},...
                                       [-0.07,0.07]);  
                                   
        [~, ~, ~, sig_t_maxClus_Abs_var]    = sign_flip_AcrossTime_PS(Abs_var,10000,[],[100,0],true,0.05,false); 
                                   
        [sig_where] =  mk_clusterStats(Abs_var,...
                                       sig_t_maxClus_Abs_var,...
                                       true,false,size(Stable_var,1)-1,0.05,true,false,[9],'Individual Seqs',...
                                       {'Absent - present'},...
                                       [-0.07,0.07]); 
    end
      
%         save(file_name,...
%              'Distvary_BB_all','stable_BB_all','stable_BB_wrong_all','vary_BB_all','vary_BB_wrong_all',...
%              'BB_tw_all','BB_wrong_tw_all','distances_seq',...
%              'Corr_Trials',...
%              'RTs',...
%              'preds_present','preds_absent',...
%              'preds_all_allSub',...
%              '-v7.3')

end