% Obtain pairwise sequenceness for inference period
% GLM approach

function [] = Replay_InferenceTime_3Length(which_data,do_replayAna,do_localiser,...
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
%         TS_RELprevBest = 30;   
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
%     idx_incr = 25; % file_name = 'idx_incr25_Length3_Replay_InferenceTime_SepIntervals.mat';

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
        file_name = 'Length3_Replay_InferenceTime.mat'; % sliding window
%         file_name = 'Control_Length3_Replay_InferenceTime.mat'; % sliding window
    elseif idx_incr == 50
        file_name = 'Length3_Replay_InferenceTime_SepIntervals.mat'; % separate intervals
%         file_name = 'Control_Length3_Replay_InferenceTime_SepIntervals.mat'; % separate intervals
    elseif idx_incr == 25
        file_name = 'idx_incr25_Length3_Replay_InferenceTime_SepIntervals.mat'; % separate intervals
%         file_name = 'Control_idx_incr25_Length3_Replay_InferenceTime_SepIntervals.mat'; % separate intervals
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness       
    BB_tw_all       = nan(length(subject),maxLag,length(tw_all),4);
    
    if do_replayAna
    
        for idx_timeW = 1:length(tw_all)

            tw = tw_all{idx_timeW};      

            s_ElStabEl       = nan(length(subject),maxLag,trial_max,n_perm,4);

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
                       
                       nSamples = size(X,1);

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
                   
                   % Prepare T matrix    
                   [T_ElEl,~,~,T_EE_wrong,~,~,~,~,~] = mk_experienced_transitions(bricks_conn_trial(idx_trial,:),bricks_rel_trial(idx_trial,:),nstates);
                   present_BB = setdiff(el_present,4);
                   
                   
                   BB_alwaysPrsnt = find([all(sum(ismember(bricks_conn_trial,1),2)) all(sum(ismember(bricks_conn_trial,3),2)) ...
                                          all(sum(ismember(bricks_conn_trial,3),2)) all(sum(ismember(bricks_conn_trial,4),2))]); 

                   BB_vary = setdiff(1:4,BB_alwaysPrsnt);
                   T_ElEl_varyTostable    = T_ElEl; T_ElEl_varyTostable(:,BB_vary) = 0; % present -> stable
                   T_ElEl_varyTovary      = T_ElEl; 
                   T_ElEl_varyTovary(BB_alwaysPrsnt,:) = 0; 
                   T_ElEl_varyTovary(:,BB_alwaysPrsnt) = 0; % present <-> present  
                   
                   if any(T_ElEl_varyTovary(:))==1
                       distantTostable = setdiff(find(sum(T_ElEl_varyTovary,2)),find(sum(T_ElEl_varyTostable,2)));                  
                   end
                   
                   % present - stable - present
                   Tpsp_Y =[present_BB(2),present_BB(1)]; % third point?
                   Tpsp_X2=[4,4]; % second point?
                   Tpsp_X1=[present_BB(1),present_BB(2)]; % first point?

                   Tpsp2=zeros(length(Tpsp_Y),nstates);
                   Tpspauto=zeros(length(Tpsp_Y),nstates);
                   for ist=1:length(Tpsp_X1)
                       Tpsp2(ist,Tpsp_Y(ist))=1;
                       Tpspauto(ist,unique([Tpsp_X1(ist),Tpsp_X2(ist)]))=1;
                   end
                   
                   % present - present - stable
                   Tpps_Y =[4,4]; % third point?
                   Tpps_X2=[present_BB(2),present_BB(1)]; % second point?
                   Tpps_X1=[present_BB(1),present_BB(2)]; % first point?

                   Tpps2=zeros(length(Tpps_Y),nstates);
                   Tppsauto=zeros(length(Tpps_Y),nstates);
                   for ist=1:length(Tpps_X1)
                       Tpps2(ist,Tpps_Y(ist))=1;
                       Tppsauto(ist,unique([Tpps_X1(ist),Tpps_X2(ist)]))=1;
                   end
                   
                   if any(T_ElEl_varyTovary(:))==1
                       Tdps_Y =[4]; % third point?
                       Tdps_X2=[setdiff(present_BB,distantTostable)]; % second point?
                       Tdps_X1=[distantTostable]; % first point?

                       Tdps2=zeros(length(Tdps_Y),nstates);
                       Tdpsauto=zeros(length(Tdps_Y),nstates);
                       for ist=1:length(Tdps_X1)
                           Tdps2(ist,Tdps_Y(ist))=1;
                           Tdpsauto(ist,unique([Tdps_X1(ist),Tdps_X2(ist)]))=1;
                       end
                       
                       Tpds_Y =[4]; % third point?
                       Tpds_X2=[distantTostable]; % second point?
                       Tpds_X1=[setdiff(present_BB,distantTostable)]; % first point?

                       Tpds2=zeros(length(Tpds_Y),nstates);
                       Tpdsauto=zeros(length(Tpds_Y),nstates);
                       for ist=1:length(Tpds_X1)
                           Tpds2(ist,Tpds_Y(ist))=1;
                           Tpdsauto(ist,unique([Tpds_X1(ist),Tpds_X2(ist)]))=1;
                       end
                   end

                   % Core sequence detection
                   X=preds_all;
                   Y=X;

                   X2bin=nan(maxLag, nSamples, nstates, nstates);  

                   for ilag=1:maxLag
                       pad = zeros([ilag nstates]);         
                       X1 = [pad; pad; X(1:end-2*ilag,:)];
                       X2 = [pad; X(1:end-ilag,:)];

                       for i=1:nstates
                           for j=1:nstates
                               X2bin(ilag,:,i,j)=X1(:,i).*X2(:,j);
                           end
                       end
                   end

                   betaPSP=nan(maxLag,length(Tpsp_Y),nstates);
                   betaPPS=nan(maxLag,length(Tpps_Y),nstates);
                   if any(T_ElEl_varyTovary(:))==1
                       betaDPS=nan(maxLag,length(Tdps_Y),nstates);
                       betaPDS=nan(maxLag,length(Tpds_Y),nstates);
                   end
                   
                   for ilag=1:maxLag        
                       Xmatrix=squeeze(X2bin(ilag,:,:,:));

                       for istate=1:length(Tpsp_Y)
                           Xpsp=squeeze(Xmatrix(:,:,Tpsp_X2(istate)));
                           Xpps=squeeze(Xmatrix(:,:,Tpps_X2(istate)));
                
                           pad  = zeros([ilag nstates]); 
                           X_control = [pad; X(1:end-ilag,:)];                           

                           tempF = pinv([Xpsp X_control ones(length(Xpsp),1)])*Y; 
                           betaPSP(ilag,istate,:)=tempF(Tpsp_X1(istate),:); 
                           
                           tempF = pinv([Xpps X_control ones(length(Xpps),1)])*Y;                
                           betaPPS(ilag,istate,:)=tempF(Tpps_X1(istate),:); 
                           
                       end
                       
                       if any(T_ElEl_varyTovary(:))==1
                           for istate=1:length(Tpds_Y)
                               Xdps=squeeze(Xmatrix(:,:,Tdps_X2(istate)));
                               Xpds=squeeze(Xmatrix(:,:,Tpds_X2(istate)));
                               
                               pad  = zeros([ilag nstates]);
                               X_control = [pad; X(1:end-ilag,:)];                               
                               
                               tempF = pinv([Xdps X_control ones(length(Xdps),1)])*Y;
                               betaDPS(ilag,istate,:)=tempF(Tdps_X1(istate),:); 
                               
                               tempF = pinv([Xpds X_control ones(length(Xpds),1)])*Y;
                               betaPDS(ilag,istate,:)=tempF(Tpds_X1(istate),:); 
                           end
                       end
                       
                   end

                   betaPSP=reshape(betaPSP,[maxLag,length(Tpsp_Y)*nstates]);
                   betaPPS=reshape(betaPPS,[maxLag,length(Tpps_Y)*nstates]);
                   if any(T_ElEl_varyTovary(:))==1
                       betaDPS=reshape(betaDPS,[maxLag,length(Tdps_Y)*nstates]);
                       betaPDS=reshape(betaPDS,[maxLag,length(Tpds_Y)*nstates]);
                   end
                   
                   % if stable building block is in the middle
                   constF=ones(length(Tpsp_Y),nstates);     

                   cc=pinv([squash(Tpsp2) squash(Tpspauto) squash(constF)])*betaPSP'; 
                   
                   s_ElStabEl(idx_sub,:,idx_trial,1,1)    = cc(1,:);
                   
                   if any(T_ElEl_varyTovary(:))==1
                       % distant present - present - stable
                       constF=ones(length(Tdps_Y),nstates);     
                       
                       cc=pinv([squash(Tdps2) squash(Tdpsauto) squash(constF)])*betaDPS'; 
                       
                       s_ElStabEl(idx_sub,:,idx_trial,1,2)    = cc(1,:);
                       
                       % present - distant present - stable
                       constF=ones(length(Tpds_Y),nstates);     
                       
                       cc=pinv([squash(Tpds2) squash(Tpdsauto) squash(constF)])*betaPDS'; 
                       
                       s_ElStabEl(idx_sub,:,idx_trial,1,3)    = cc(1,:);
                   end
                   
                   % if stable building block is at the end
                   constF=ones(length(Tpps_Y),nstates);     

                   cc=pinv([squash(Tpps2) squash(Tppsauto) squash(constF)])*betaPPS'; 
                   
                   s_ElStabEl(idx_sub,:,idx_trial,1,4)    = cc(1,:);

                end               

            end

            fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))
            
            s_ElStabEl       = mk_MeanSess(s_ElStabEl);
            
            s_ElStabEl       = mk_CombineSess(s_ElStabEl,combine_sess); 
            
            fprintf('########Done with tw %d of %d.########\n',idx_timeW,length(tw_all))
        
        end
    
    else

        load(file_name)

    end         

    var_stab_var   = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,1),4),2)); % absent -> vary
    dvar_var_stab  = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,2),4),2)); % absent -> vary
    var_dvar_stab  = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,3),4),2)); % absent -> vary
    var_var_stab   = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,4),4),2)); % absent -> vary      

    if size(var_stab_var,2) <= 50

        [h_VarStabVar, ~]  = ttest(var_stab_var);
        [h_DVarVarStab, ~] = ttest(dvar_var_stab);
        [h_VarDVarStab, ~] = ttest(var_dvar_stab);
        [h_VarVarStab, ~]  = ttest(var_var_stab);
    
        if size(var_stab_var,2) == 7
            mk_PlotSotErrbar_Group_PS(cat(3,var_stab_var,dvar_var_stab,var_dvar_stab,var_var_stab),true,'Individual Sequences',...
                                      'Time-Interval','Length 3 Different Time-Intervals',...
                                      {'-500ms - 500ms' '0ms - 1000ms' '500ms - 1500ms' '1000ms - 2000ms' '1500ms - 2500ms' '2000ms - 2500ms' '2500ms - 3500ms'},...
                                      cat(3,h_VarStabVar,h_DVarVarStab,h_VarDVarStab,h_VarVarStab),[1,2,3,4],...
                                      {'[present -> stable] -> present' '[distant present -> present] -> stable' '[present -> distant present] -> stable' '[present -> present] -> stable'},...
                                      [-0.5,0.5],[0 8],[])
        elseif size(var_stab_var,2) == 13
            mk_PlotSotErrbar_Group_PS(cat(3,var_stab_var,dvar_var_stab,var_dvar_stab,var_var_stab),true,'Individual Sequences',...
                                      'Time-Interval','Length 3 Different Time-Intervals',...
                                      {'-500ms - 500ms' '-250ms - 750ms' '0ms - 1000ms' '250ms - 1250ms' '500ms - 1500ms' '750ms - 1750ms' '1000ms - 2000ms' '1250ms - 2250ms' '1500ms - 2500ms' '1750ms - 2750ms'  '2000ms - 3000ms'  '2250ms - 3250ms'  '2500ms - 3500ms'},...
                                      cat(3,h_VarStabVar,h_DVarVarStab,h_VarDVarStab,h_VarVarStab),[1,2,3,4],...
                                      {'[present -> stable] -> present' '[distant present -> present] -> stable' '[present -> distant present] -> stable' '[present -> present] -> stable'},...
                                      [-0.5,0.5],[0 14],[])
        end                            
                              
    else
        % plot individual variables
        [~, ~, ~, sig_t_maxClus_VarStabVar]     = sign_flip_AcrossTime_PS(var_stab_var,10000,[],[100,0],true,0.05,false); 
        [~, ~, ~, sig_t_maxClus_dvar_var_stab]  = sign_flip_AcrossTime_PS(dvar_var_stab,10000,[],[100,0],true,0.05,false); 
        [~, ~, ~, sig_t_maxClus_var_dvar_stab]  = sign_flip_AcrossTime_PS(var_dvar_stab,10000,[],[100,0],true,0.05,false);
        [~, ~, ~, sig_t_maxClus_var_var_stab]   = sign_flip_AcrossTime_PS(var_var_stab,10000,[],[100,0],true,0.05,false);

        [sig_where] =  mk_clusterStats(cat(3,var_stab_var,dvar_var_stab,var_dvar_stab,var_var_stab),...
                                       cat(3,sig_t_maxClus_VarStabVar,sig_t_maxClus_dvar_var_stab,sig_t_maxClus_var_dvar_stab,sig_t_maxClus_var_var_stab),...
                                       true,false,size(var_stab_var,1)-1,0.05,true,false,[1,2,3,4],'Individual Sequences',...
                                       {'[present -> stable] -> present' '[distant present -> present] -> stable' '[present -> distant present] -> stable' '[present -> present] -> stable'},...
                                       []); 
        time = -50:350;
        disp(sig_where)
        time(sig_where{1}{2})
        time(sig_where{1}{3})
        time(sig_where{1}{4})
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