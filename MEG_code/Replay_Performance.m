% Obtain pairwise sequenceness for inference period
% GLM approach

function [] = Replay_Performance()
    %% Prelim

    idx_incr = 1; % this reproduces main finding in figure 6
%     idx_incr = 50; % this reproduces supplementary figure 5


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% some relevant variables    
    % Filename of pre-analysed sequence data
    if idx_incr == 1
        file_name = 'Replay_InferenceTime_Performance.mat'; % sliding window
    elseif idx_incr == 50
        file_name = 'Replay_InferenceTime_SepIntervals_Performance.mat'; % separate intervals
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain different sequenceness    
    load(file_name)       

    Abs_var    = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,1),4),2)); % absent -> vary
    Var_abs    = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,2),4),2)); % vary -> absent
    Abs_stable = squeeze(nanmean(nanmean(BB_wrong_tw_all(:,:,:,3),4),2)); % absent -> stable 
    
    Var_var    = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,1),4),2)); % vary -> vary
    Var_stable = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,2),4),2)); % vary -> stable
    Stable_var = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,3),4),2)); % stable -> vary
    DistVar_stable = squeeze(nanmean(nanmean(BB_tw_all(:,:,:,4),4),2)); % stable -> vary        
    
    mk_PlotSotErrbar_Group_PS(RTs')
    mk_PlotSotErrbar_Group_PS(Corr_Trials')
    
    RT_mean = nanmean(RTs,2);
    Corr_mean = nanmean(Corr_Trials,2);
    
    %% effects for mean across subjects
    
    
    % sig individual effects                         
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
   
   % first peak only
   sig_timeWin_mean = [nanmean(Abs_stable(:,sig_where{1}{3}),2), ...
                       nanmean(DistVar_stable(:,sig_where{1}{4}),2), ...
                       nanmean(Var_stable(:,sig_where{1}{5}(1:45)),2)];
                   
   [r,p] = corr([mean(sig_timeWin_mean,2),RT_mean,Corr_mean])
   
   % last emergence of present-stable
   sig_timeWin_mean = [nanmean(Var_stable(:,sig_where{1}{5}(46:end)),2)];
   
   [r,p] = corr([sig_timeWin_mean,RT_mean,Corr_mean])
   
   %% trial-by-trial effects
     
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

   beta_subs_all_RT   = zeros(size(s_ElEl_allTrls,1),2);
   beta_subs_all_Corr = zeros(size(s_ElEl_allTrls,1),2);
                                   
   for idx_sub = 1:size(s_ElEl_allTrls,1)
       
       RTs_sub         = RTs(idx_sub,:);
       Corr_Trials_sub = Corr_Trials(idx_sub,:);
       
       idx_use = ~isnan(RTs_sub);
       
       RTs_sub         = RTs_sub(idx_use);
       Corr_Trials_sub = Corr_Trials(idx_use);
       
       Abs_stable     = squeeze(s_ElEl_wrong_allTrls(idx_sub,:,3,:));
       DistVar_stable = squeeze(s_ElEl_allTrls(idx_sub,:,4,:));
       Var_stable     = squeeze(s_ElEl_allTrls(idx_sub,:,2,:));
       
       Abs_stable       = Abs_stable(:,sig_where{1}{3});
       DistVar_stable   = DistVar_stable(:,sig_where{1}{4});
       Var_stable_1Peak = Var_stable(:,sig_where{1}{5}(1:45));
       Var_stable_2Peak = Var_stable(:,sig_where{1}{5}(46:end));

       dm = [nanmean([nanmean(Abs_stable,2) nanmean(DistVar_stable,2) nanmean(Var_stable_1Peak,2)],2) ...
             nanmean(Var_stable_2Peak,2)];
         
       dm = dm(idx_use,:);
       
       use_trls = (isnan(dm(:,1))==0 & isnan(dm(:,2))==0);
       
       RTs_sub         = RTs_sub(use_trls);
       Corr_Trials_sub = Corr_Trials_sub(use_trls);
       
       dm = dm(use_trls,:);
       
       cm = eye(size(dm,2)+1);
       
       [beta,~,~] = ols_PS(RTs_sub',[dm ones(size(dm,1),1)],cm);       
       
       beta_subs_all_RT(idx_sub,:) = beta(1:end-1);
       
       [b,dev,stats] = glmfit(dm,Corr_Trials_sub','binomial');
       
       beta_subs_all_Corr(idx_sub,:) = b(2:end);
       
   end
   
    
end