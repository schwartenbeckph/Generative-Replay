% Obtain pairwise sequenceness for inference period
% GLM approach

function [] = Replay_PlotClassifier(which_data,do_replayAna,do_localiser,...
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

    idx_tw = 1;
    tw_all{idx_tw} = start_tw:end_tw;
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
    
    idx_trialbeg = abs(trlbeg)+1;
    
    trial_max = 6*48;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    file_name = 'Class_preds_inference';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now obtain classifier predictions  
    
    if do_replayAna
    
        for idx_timeW = 1:length(tw_all)

            tw = tw_all{idx_timeW};      

            Corr_Trials = nan(length(subject),trial_max);
            RTs         = nan(length(subject),trial_max);        

            preds_all_allSub = nan(length(subject),length(tw),nstates);

            for idx_sub = 1:length(subject)   

                behav_folder   = fullfile(behav_result_path,(subject{idx_sub}));

                %%%%% get biased brick info %%%%%
                [~,~,~,bricks_rel_trial,~,~,correct_trials_all,rt_all,~,~] = mk_bricks_sub(behav_folder); % obtain brick info of subject based on observed stims

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
    %             figure,plot(ERF_subj(idx_sub,:))

                for idx_trial = use_trials   

                   if isempty(tw)
                       X = data(:,idx_trialbeg:end,idx_trial);
                   elseif length(tw)==1
                       X = data(:,idx_trialbeg+tw:end,idx_trial);
                   else
                       X = data(:,idx_trialbeg+tw,idx_trial);
                   end
                   
%                    X = data(:,:,idx_trial);

                   X = nanmean(X,3)';            

                   if ~exist('preds_present')
                       preds_present = nan(size(X,1),length(subject),trial_max);
                       preds_absent = nan(size(X,1),length(subject),trial_max);
                       preds_stable = nan(size(X,1),length(subject),trial_max);
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

                       preds_present(:,idx_sub,idx_trial) = mean(preds_el(:,setdiff(el_present,4)),2);
                       preds_absent(:,idx_sub,idx_trial)  = mean(preds_el(:,el_absent),2);
                       preds_stable(:,idx_sub,idx_trial)  = mean(preds_el(:,4),2);

                   end

                end  
                
                fprintf('####Subject %d of %d done.####\n',idx_sub,length(subject))

            end

        end
    
    else

        load(file_name)

    end      
    
    %% now plot    
    col   = {[0, 0.4470, 0.7410], ...       % blue
             [0.4940, 0.1840, 0.5560], ...  % purple
             [0.4660, 0.6740, 0.1880], ...  % green
             [0.6350, 0.0780, 0.1840], ...  % red
             [0.3010, 0.7450, 0.9330], ...  % cyan
             [0.9290, 0.6940, 0.1250], ...  % yellow
             [0.0000, 0.2600, 0.1500], ...  % 'british racing green'
             [0.8500, 0.3250, 0.0980], ...  % orange
             [0, 0, 0], ...                 % black
             [0.6, 0.6, 0.6]};              % grey  
    
    % Now: plt present, absent, stable
    greys_col = {[0.3, 0.3, 0.3],...
                 [0.5, 0.5, 0.5],...
                 [0.7, 0.7, 0.7]};
             
    preds_all_allSub_plot = cat(3,nanmean(preds_stable,3)',nanmean(preds_absent,3)',nanmean(preds_present,3)');
         
    n_time = size(preds_all_allSub_plot,2);
    n_bricks = size(preds_all_allSub_plot,3);
    
    figure, set(gcf,'color','white'),hold on
    
    col_idx = [4 5 3];
    
    for idx_data = [1:n_bricks]
    
        data_plot = preds_all_allSub_plot(:,:,idx_data);

        mm = nanmean(data_plot,1);
        ss = nanstd(data_plot,0,1)./sqrt(size(data_plot,1));

        p(idx_data) = plot(mm,'Color', col{col_idx(idx_data)},'LineWidth',3,'LineStyle',':');
        fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{col_idx(idx_data)},'EdgeAlpha',0,'FaceAlpha',0.2);
    
    end
    
    % plot zone where candidate->stable replay effects become
    % significant
    xline(idx_trialbeg,'--','Color',col{end-1})
    xline(idx_trialbeg+18,'--','Color',col{end-2})
    xline(idx_trialbeg+69,'--','Color',col{end-2})
    a = area([idx_trialbeg+18,idx_trialbeg+69],[0.4,0.4],'FaceAlpha',0.2);
    a(1).FaceColor = col{end-2};
    a(1).EdgeColor = col{end-2};
    
    % plot zone where present<->absent vs. present<->present replay effects become
    % significant
    xline(idx_trialbeg+27,'--','Color',col{end-3})
    xline(idx_trialbeg+66,'--','Color',col{end-3})
    a = area([idx_trialbeg+27,idx_trialbeg+66],[0.4,0.4],'FaceAlpha',0.2);
    a(1).FaceColor = col{end-3};
    a(1).EdgeColor = col{end-3};  
    
    % plot zone where present->stable replay effects become
    % significant
    xline(idx_trialbeg+160,'--','Color',col{3})
    a = area([idx_trialbeg+160,n_time],[0.4,0.4],'FaceAlpha',0.2);
    a(1).FaceColor = col{3};
    a(1).EdgeColor = col{3};    
    
    legend(p,"Stable","Absent","Candidate")
    
    ylabel('')
    xlabel('')
    xticks([idx_trialbeg idx_trialbeg+100 idx_trialbeg+200 idx_trialbeg+300])
    xticklabels({'' '' '' ''})
    yticks([0 0.1 0.2 0.3])
    yticklabels({'' '' '' ''})
    
    xlim([0,401])
    ylim([0,0.3])   
    hold off

end