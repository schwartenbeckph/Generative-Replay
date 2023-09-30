% This script obtains individual accuracies from the element/relation
% decoding under different levels of smoothing and does group stats on
% them.

function [overall_accuracy, overall_accuracy_bestTS, se_accuracy_bestTS, L1_best_subj, TS_best_subj, overall_accuracy_L1] = ...
    AccuracyDecoding_allTS(scan_result_path,decode,do_plot)

    %% prelim
    if nargin==0

        based = '/Users/Philipp/Dropbox/StimConstrMEG/Results_Scanning/code_MEG_Final';

        scan_result_path    = fullfile(based,'data');

        decode = 'localiser';

        do_plot = true;

    end

    subject = dir(fullfile(scan_result_path,'s*'));
    subject = {subject.name};

    nstates = 4; % always four classes - either elements or relations

    chance_level = 1/nstates;

    ENL1          = 0.001:0.001:0.01;
    
    ENL1_idx_use = 1:length(ENL1);
    ENL1 = ENL1(ENL1_idx_use);
    %% get accuracies
    L1_best = zeros(1,length(ENL1));

    L1_best_subj        = nan(length(subject),1);
    TS_best_subj        = nan(length(subject),1);
    MaxVal_TS_best_subj = nan(length(subject),1);

    idx_data = 1;

    %%%%% get data %%%%%
    file_name = ['DataDecode_ALLTS_',decode];
    
    load(file_name) % load 'accur' and 'preds'
    
    n_sess = size(accur_all,5);
    
    for idx_sess=1:n_sess

        for idx_sub = 1:length(subject)

            %%%%%  obtain accuracies, predictions %%%%% 
            accur = squeeze(accur_all(idx_sub,ENL1_idx_use,:,:,idx_sess));

            [~,L1_max] = max(max(max(accur,[],3),[],2));

            accuracy_all(idx_data,:,:)  = accur(L1_max,:,:);

            [max_decoding,max_timepoint] = max(max(accur(L1_max,:,:)));
            accuracy_all_bestTS(idx_data,:) = accur(L1_max,:,max_timepoint);

            L1_best(L1_max) = L1_best(L1_max)+1;

            L1_best_subj(idx_sub) = ENL1(L1_max);

            TS_best_subj(idx_sub)        = max_timepoint;
            MaxVal_TS_best_subj(idx_sub) = max_decoding;

            idx_data = idx_data + 1;   

            fprintf('Subject %d of %d done.\n',idx_sub,length(subject))

        end

        [~,overall_accuracy_L1] = max(max(max(accur,[],3),[],2));

        overall_accuracy_L1 = ENL1(overall_accuracy_L1);

        overall_accuracy        = squeeze(nanmean(accuracy_all,1));

        [~, max_timepoint_Mean] = max(max(overall_accuracy,[],2));

        max_data = squeeze(accuracy_all(:,max_timepoint_Mean,:));
        
        overall_accuracy_bestTS = nanmean(max_data,1);
        se_accuracy_bestTS      = nanstd(max_data,0,1)./sqrt(size(max_data,1));

        %% plot

        if do_plot

            tss = 0:length(overall_accuracy_bestTS)-1;

            col   = {[0, 0.4470, 0.7410], ...       % blue
                     [0.4940, 0.1840, 0.5560], ...  % purple
                     [0.4660, 0.6740, 0.1880], ...  % green
                     [0.9350, 0.1780, 0.2840], ...  % red
                     [0.3010, 0.7450, 0.9330], ...  % cyan
                     [0.9290, 0.6940, 0.1250], ...  % yellow
                     [0, 0, 0], ...                 % black
                     [0.6, 0.6, 0.6]};              % grey


            %%%%% plot all decoding %%%%%
            figure
            set(gcf,'color','white')
            imagesc(overall_accuracy),colorbar,caxis([0 max(max(overall_accuracy))])

            xlabel('Time in 10ms after onset (tested)'),ylabel('Time in 10ms after onset (trained)')  

            %%%%% plot all decoding best TS %%%%%
            figure,hold on
            set(gcf,'color','white')
            
            overall_accuracy_bestTS = nanmean(accuracy_all(:,:,max_timepoint_Mean),1);

            plot(tss,overall_accuracy_bestTS,'Color',col{7},'LineWidth',3);     
            
            fill([tss, fliplr(tss)],[overall_accuracy_bestTS+se_accuracy_bestTS, fliplr(overall_accuracy_bestTS-se_accuracy_bestTS)], col{8},'EdgeAlpha',0,'FaceAlpha',0.2);
            
            line([max_timepoint_Mean-1 max_timepoint_Mean-1], [0 1], 'Color',[.8 .8 .8],'LineWidth',3); % -1 because we start at 0!  

            xlabel('Time in 10ms after onset'),ylabel('Accuracy')
            set(gca, 'XTick', [tss(1):10:tss(end)])
            set(gca, 'XTickLabel', tss(1:10:end))
            ylim([0,0.7])
            hold off

            %%%%% plot individual decoding best TS %%%%%
            if idx_sess==1 || idx_sess==n_sess
                figure
                set(gcf,'color','white')

                for idx_sub = 1:length(subject)

                    subplot(4,5,idx_sub),hold on
                    plot(tss,accuracy_all_bestTS(idx_sub,:))

                    xlabel('Time in 10ms before/after onset'),ylabel('Accuracy')
                    
                    title(sprintf('Individual decoding sess %d',idx_sess))
                    
                    plot(tss,ones(length(accur),1)*chance_level,'--k','MarkerSize',12)
                    axis([tss(1) tss(end) 0 1])
                    set(gca, 'XTick', tss(1):20:tss(end))
                    hold off

                end    
            end

            %%%%% plot best L1 %%%%%
            if idx_sess==1 || idx_sess==n_sess
                figure
                set(gcf,'color','white')
                bar(L1_best)
                xlabel('L1 penalty'),ylabel('Frequency')
                set(gca, 'XTick', 1:length(ENL1))
                set(gca, 'XTickLabel', ENL1)
                ylim([0,length(subject)])
                set(gca, 'YTick', 0:2:length(subject))
                title(sprintf('Best L1s sess %d',idx_sess))
            end

            %%%%% plot best TS %%%%%
            if idx_sess==1 || idx_sess==n_sess
                figure
                set(gcf,'color','white')
                hist(tss(TS_best_subj))
                xlabel('Time step training'),ylabel('Frequency')
                axis([tss(1) tss(end) 0 10])
                set(gca, 'XTick', tss(1):10:tss(end))
                title(sprintf('Best training time step for decoding sess %d',idx_sess)) 
            end

        end
    
    end

end


