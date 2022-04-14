% This function plots the sequenceness across sessions and differences
% between those sessions
% 
% if two sequence data are entered this is done for the difference between
% those measures for every session
% 
% S is #subjects x time points x sessions
function [] = plot_diff_sequenceness(S1,S2,title_plot,n_perm,do_reverse,do_one_sided,do_shuffle,legend_lab,type,plot_se,tell_sig,covariate,title_covariate,...
                                     diff_plot,no_shuffrange,which_data,no_stats,col_idx)

    %% prelimary stuff                  
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
         
%              [0.9350, 0.1780, 0.2840], ...  % red          
    
    if nargin<2
        S2           = [];
        S_diff       = S1;
        n_perm       = 10000;
        title_plot   = [];
        do_reverse   = true;
        do_one_sided = false;
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<3
        n_perm       = 10000;
        title_plot   = [];  
        do_reverse   = true;
        do_one_sided = false;
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<4
        n_perm       = 10000;
        do_reverse   = true;
        do_one_sided = false;
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<5
        do_reverse   = true;
        do_one_sided = false;
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<6
        do_one_sided = false;
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<7
        do_shuffle   = false;
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<8
        legend_lab   = [];
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<9
        type         = 'mean';
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<10
        plot_se      = true;
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<11
        tell_sig     = false;
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<12
        covariate    = [];
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<13
        title_covariate = [];
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<14
        diff_plot    = false;
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<15
        no_shuffrange = false;
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<16
        which_data   = 'plan';
        no_stats     = false;
        col_idx      = [];
    elseif nargin<17
        no_stats     = false;
        col_idx      = [];
    elseif nargin<18
        col_idx      = [];
    end
    
    if isempty(S2)
        S_diff = S1;
    else
        if size(S1,1)~=size(S2,1) || size(S1,2)~=size(S2,2) || size(S1,3)~=size(S2,3)
            error('I cannot deal with this')
        end
        S_diff = S1 - S2;
    end    
    
    n_subs = size(S1,1);
    n_time = size(S1,2);
    n_sess = size(S1,3);
    n_shuf = size(S1,4);
    n_data = size(S1,5);
    
    if isempty(col_idx)
       col_idx = 1:n_data; 
    end
    
    if do_one_sided
%         perc_bound_up   = 95;
        perc_bound_up    = 100;
        perc_bound_down = 0;
    else
%         perc_bound_up   = 97.5;
        perc_bound_up   = 100;
%         perc_bound_down = 2.5;
        perc_bound_down = 0;
    end
    
    for idx_data=1:n_data
    
        %% obtain permutation stats per session
%         if do_shuffle % || ~isempty(S2)
        if n_shuf>1 || ~do_shuffle

            up_sess  = zeros(n_sess,n_time);
            low_sess = zeros(n_sess,n_time);
            emp_sess = zeros(n_sess,n_time);

            if ~do_shuffle

                for idx_sess=1:n_sess

                    [up_sess(idx_sess,:), low_sess(idx_sess,:), emp_sess(idx_sess,:)] = sign_flip_AcrossTime_PS(S_diff(:,:,idx_sess,idx_data),n_perm);

                end

            else

                for idx_sess=1:n_sess

                    [up_sess(idx_sess,:), low_sess(idx_sess,:), emp_sess(idx_sess,:)] = sign_flip_AcrossTime_PS([],[],squeeze(S_diff(:,:,idx_sess,:,idx_data)));

                end

            end

        else

            up_sess  = [];
            low_sess = [];
            emp_sess = [];

        end

        %% obtain permutation stats for differences between sessions
        if n_sess==1
            all_pairwise_diff = [1 1];
        else
            if do_reverse
                all_pairwise_diff = fliplr(nchoosek(1:n_sess,2));
            else
                all_pairwise_diff = nchoosek(1:n_sess,2);
            end
        end

        if size(all_pairwise_diff,1)>8
            all_pairwise_diff = all_pairwise_diff(1:2:end,:);
        end

        up_sess_diff  = zeros(size(all_pairwise_diff,1),n_time);
        low_sess_diff = zeros(size(all_pairwise_diff,1),n_time);
        emp_sess_diff = zeros(size(all_pairwise_diff,1),n_time);

        S_diff_sess = zeros(n_subs,n_time,size(all_pairwise_diff,1),n_shuf,n_data);

        for idx_sess=1:size(all_pairwise_diff,1)

            S_diff_sess(:,:,idx_sess,:,:) = S_diff(:,:,all_pairwise_diff(idx_sess,1),:,:) - S_diff(:,:,all_pairwise_diff(idx_sess,2),:,:);

            if ~do_shuffle

                [up_sess_diff(idx_sess,:), low_sess_diff(idx_sess,:), emp_sess_diff(idx_sess,:)] = sign_flip_AcrossTime_PS(S_diff_sess(:,:,idx_sess,idx_data),n_perm);

            else

                [up_sess_diff(idx_sess,:), low_sess_diff(idx_sess,:), emp_sess_diff(idx_sess,:)] = sign_flip_AcrossTime_PS([],[],S_diff_sess(:,:,idx_sess,:,idx_data));

            end

        end

        %% now plot
        S_diff_plot      = S_diff(:,:,:,1,idx_data);
        S_diff_sess_plot = S_diff_sess(:,:,:,1,idx_data); 

        % Figure 1: all individual sessions
        if strcmp(type,'mean')
            
            if ~isempty(up_sess)
                max_min_data = max([max(max(squeeze(nanmean(S_diff_plot,1)))), abs(min(min(squeeze(nanmean(S_diff_plot,1)))))]);

                y_max = max([max(max(up_sess)) abs(min(min(low_sess)))]);

                y_max = max(y_max,max_min_data);

                y_max = round(y_max*2,2);
            else
                y_max = max(max(squeeze(nanmean(S_diff_plot,1))));
                y_max = round(y_max*3,2);
            end
            
            if idx_data==1
                y_max_all = y_max;
            else
                y_max_all = max(y_max_all,y_max);
            end

            if idx_data==1
                figure, set(gcf,'color','white')
            end
            if n_sess>3
                n_rows = 2;
            else
                n_rows = 1;
            end
            n_cols = n_sess;

            for idx_sess=1:size(S_diff_plot,3)
                subplot(n_rows,n_cols,idx_sess),hold on

                data_sess = S_diff_plot(:,:,idx_sess);

                mm = nanmean(data_sess,1);
                ss = nanstd(data_sess,0,1)./sqrt(size(data_sess,1));

                p(idx_data) = plot(mm,'Color',col{col_idx(idx_data)});
                if plot_se
                    fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{col_idx(idx_data)},'EdgeAlpha',0,'FaceAlpha',0.2);
                end

%                 if y_max~=0
%                     ylim([-y_max,y_max])
%                 end
                
                if y_max_all~=0
                    ylim([-y_max_all,y_max_all])
                end

                if ~isempty(up_sess) && ~no_stats
                    if n_data==1 && ~no_shuffrange
                        fill([1:n_time, fliplr(1:n_time)],[up_sess(idx_sess,:), fliplr(low_sess(idx_sess,:))], col{8},'EdgeAlpha',0,'FaceAlpha',0.2);
                    end
                    if do_one_sided
                        if n_data==1
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(idx_sess,:),perc_bound_up),'--k','MarkerSize',6) % one sided test
                        else
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(idx_sess,:),perc_bound_up),'--','Color',col{col_idx(idx_data)},'MarkerSize',6) % one sided test
                        end
                        if tell_sig
                            sig_idx = find(mm>=prctile(up_sess(idx_sess,:),perc_bound_up));
                            if ~isempty(sig_idx)
                                disp(['Sig at ',num2str(sig_idx)])
                                disp(['Max at ',num2str(find(mm==max(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                            end
                        end
                    else
                        if n_data==1
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(idx_sess,:),perc_bound_up),'--k','MarkerSize',6) % two sided test
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess(idx_sess,:),perc_bound_down),'--k','MarkerSize',6)
                        else
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(idx_sess,:),perc_bound_up),'--','Color',col{col_idx(idx_data)},'MarkerSize',6) % two sided test
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess(idx_sess,:),perc_bound_down),'--','Color',col{col_idx(idx_data)},'MarkerSize',6)
                        end
                        if tell_sig
                            sig_idx = find(mm>=prctile(up_sess(idx_sess,:),perc_bound_up));
                            if ~isempty(sig_idx)
                                disp(['Sig pos at ',num2str(sig_idx)])
                                disp(['Max at ',num2str(find(mm==max(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                                disp(['Max at ',num2str(find(mm==max(mm))),' value ',num2str(max(mm))])
                                disp(['Sig at ',num2str(prctile(up_sess(idx_sess,:),perc_bound_up))])
                            end
                            
                            sig_idx = find(mm<=prctile(low_sess(idx_sess,:),perc_bound_down));
                            if ~isempty(sig_idx)
                                disp(['Sig neg at ',num2str(sig_idx)])
                                disp(['Min at ',num2str(find(mm==min(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                                disp(['Min at ',num2str(find(mm==min(mm))),' value ',num2str(min(mm))])
                                disp(['Sig at ',num2str(prctile(low_sess(idx_sess,:),perc_bound_down))])
                            end
                        end
                    end
                end
                
                title({title_plot})
                
                if diff_plot
                    ylabel('difference')                    
                    if ~strcmp(which_data,'rest')
                        xticks([1 6:5:n_time])
                %         xticklabels({(0:5:n_time)*100})
                        xticklabels({(-5:5:n_time)*100})
                        xlabel('time interval, 500ms  starting at (ms)')
                    else
%                         xticks([1 11:10:n_time])
                        xticks(linspace(1,n_time,9))
                        xticklabels({(0:0.5:4)})
                        xlabel('time interval, 0.5 minutes  starting at (minutes)')
                    end
                else
                    ylabel('sequenceness')
                    xlabel('lag (seconds)')
                    xticks([0:10:n_time])
                    xticklabels({(0:10:n_time)/100})
                end

%                 if ~isempty(legend_lab)
%                     legend(p,legend_lab)
%                 end

                plot(1:n_time,zeros(length(1:n_time),1),'--k','MarkerSize',6)

%                 title({title_plot; ['Session ',num2str(idx_sess)]})
%                 title({title_plot})
%                 ylabel('sequenceness')
%                 xlabel('lag (seconds)')
%                 xticks([0:10:n_time])
%                 xticklabels({(0:10:n_time)/100})
            end

        % Figure 2: all session differences
        elseif strcmp(type,'diff')
            if idx_data==1
                if ~isempty(up_sess_diff)
                    y_max = max([max(max(up_sess_diff)) abs(min(min(low_sess_diff)))]);
                    y_max = round(y_max*2,2);
                else
                    y_max = max(max(squeeze(nanmean(S_diff_sess_plot,1))));
                    y_max = round(y_max*3,2);
                end
            end

            if idx_data==1
                figure, set(gcf,'color','white')
            end
            
            if size(S_diff_sess_plot,3)<=3 % rest
                n_rows = 1;
                n_cols = size(S_diff_sess_plot,3);
            else
                n_rows = 2;
                n_cols = 4;
            end

            for idx_sess=1:size(S_diff_sess_plot,3)
                subplot(n_rows,n_cols,idx_sess),hold on

                data_sess = S_diff_sess_plot(:,:,idx_sess);

                mm = nanmean(data_sess,1);
                ss = nanstd(data_sess,0,1)./sqrt(size(data_sess,1));

                p(idx_data) = plot(mm,'Color',col{col_idx(idx_data)});
                if plot_se
                    fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{col_idx(idx_data)},'EdgeAlpha',0,'FaceAlpha',0.2);
                end

                if y_max~=0
                    ylim([-y_max,y_max])
                end

                if ~isempty(up_sess) && ~no_stats
                    if n_data==1  && ~no_shuffrange 
                        fill([1:n_time, fliplr(1:n_time)],[up_sess_diff(idx_sess,:), fliplr(low_sess_diff(idx_sess,:))], col{8},'EdgeAlpha',0,'FaceAlpha',0.2);
                    end
                    if do_one_sided
                        if n_data==1
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess_diff(idx_sess,:),perc_bound_up),'--k','MarkerSize',6) % one sided test
                        else
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess_diff(idx_sess,:),perc_bound_up),'--','Color',col{col_idx(idx_data)},'MarkerSize',6) % one sided test
                        end   
                        if tell_sig
                            sig_idx = find(mm>=prctile(up_sess(idx_sess,:),perc_bound_up));
                            if ~isempty(sig_idx)
                                disp(['Sig at ',num2str(sig_idx)])
                                disp(['Max at ',num2str(find(mm==max(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                            end
                        end
                    else
                        if n_data==1
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess_diff(idx_sess,:),perc_bound_up),'--k','MarkerSize',6) % one sided test
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess_diff(idx_sess,:),perc_bound_down),'--k','MarkerSize',6)
                        else
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess_diff(idx_sess,:),perc_bound_up),'--','Color',col{col_idx(idx_data)},'MarkerSize',6) % one sided test
                            plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess_diff(idx_sess,:),perc_bound_down),'--','Color',col{col_idx(idx_data)},'MarkerSize',6)
                        end   
                        if tell_sig
                            sig_idx = find(mm>=prctile(up_sess(idx_sess,:),perc_bound_up));

                            if ~isempty(sig_idx)
                                disp(['Sig pos at ',num2str(sig_idx)])
                                disp(['Max at ',num2str(find(mm==max(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                                disp(['Max at ',num2str(find(mm==max(mm))),' value ',num2str(max(mm))])
                                disp(['Sig at ',num2str(prctile(up_sess(idx_sess,:),perc_bound_up))])
                            end

                            sig_idx = find(mm<=prctile(up_sess(idx_sess,:),perc_bound_down));                        
                            if ~isempty(sig_idx)
                                disp(['Sig neg at ',num2str(sig_idx)])
                                disp(['Min at ',num2str(find(mm==min(mm(sig_idx))))])
                            else
                                fprintf('Nothing sig.\n')
                                disp(['Min at ',num2str(find(mm==min(mm))),' value ',num2str(min(mm))])
                                disp(['Sig at ',num2str(prctile(low_sess(idx_sess,:),perc_bound_down))])
                            end
                        end                    
                    end
                end

                title({title_plot; ['Rest ',num2str(all_pairwise_diff(idx_sess,1)),' - Rest ',num2str(all_pairwise_diff(idx_sess,2))]})
%                 ylabel('sequenceness')
%                 xlabel('lag (seconds)')
%                 xticks([0:10:n_time])
%                 xticklabels({(0:10:n_time)/100})

%                 title({title_plot})
                
                if diff_plot
                    ylabel('difference')                    
                    if ~strcmp(which_data,'rest')
                        xticks([1 6:5:n_time])
                %         xticklabels({(0:5:n_time)*100})
                        xticklabels({(-5:5:n_time)*100})
                        xlabel('time interval, 1000ms  starting at (ms)')
                    else
%                         xticks([1 11:10:n_time])
                        xticks(linspace(1,n_time,9))
                        xticklabels({(0:0.5:4)})
                        xlabel('time interval, 0.5 minutes  starting at (minutes)')
                    end
                else
                    ylabel('sequenceness')
                    xlabel('lag (seconds)')
                    xticks([0:10:n_time])
                    xticklabels({(0:10:n_time)/100})
                end

%                 if ~isempty(legend_lab)
%                     legend(p,legend_lab)
%                 end
% 
                plot(1:n_time,zeros(length(1:n_time),1),'--k','MarkerSize',6)
                
                
            end
            
%             if ~isempty(legend_lab)
%                 legend(p,legend_lab)
%             end
% 
%             plot(1:n_time,zeros(length(1:n_time),1),'--k','MarkerSize',6)
        
        end
        
        
        
    
    end
    
    if ~isempty(legend_lab)
        legend(p,legend_lab)
    end

%     plot(1:n_time,zeros(length(1:n_time),1),'--k','MarkerSize',6)
    
%     plot(1:n_time,zeros(length(1:n_time),1),'--k','MarkerSize',6)
    
%     title({title_plot})
%     if diff_plot
%         ylabel('difference')
%         xlabel('time interval, 500ms  starting at (ms)')
%         xticks([1 6:5:n_time])
% %         xticklabels({(0:5:n_time)*100})
%         xticklabels({(-5:5:n_time)*100})
%     else
%         ylabel('sequenceness')
%         xlabel('lag (seconds)')
%         xticks([0:10:n_time])
%         xticklabels({(0:10:n_time)/100})
%     end
%     
%     if ~isempty(legend_lab)
%         legend(p,legend_lab)
%     end
    
    %% test for effect of covariate
    if ~isempty(covariate)
        
        [beta, var_betas,t_vals] = ols_PS(S_diff_plot,[covariate,ones(size(covariate,1),1)],eye(size(covariate,2)+1));
        
        [~, ~, ~, ~, S_diff_plot_perm] = sign_flip_AcrossTime_PS(S_diff_plot,10000);

        
        beta_perm      = nan(size(S_diff_plot_perm,3),size(beta,2));
        var_betas_perm = nan(size(S_diff_plot_perm,3),size(beta,2));
        crit_t   = tinv(.975,size(S_diff_plot,1));
        sig_t_maxClus  = nan(1,size(S_diff_plot_perm,3));
    
        for idx_perm=1:size(S_diff_plot_perm,3)
            
            [beta_temp, var_temp, t_vals_perm] = ols_PS(S_diff_plot_perm(:,:,idx_perm),[covariate,ones(size(covariate,1),1)],eye(size(covariate,2)+1));
            beta_perm(idx_perm,:) = beta_temp(1,:);
            var_betas_perm(idx_perm,:) = var_temp(1,:);
            
            t_vals_perm = t_vals_perm(1,:);
            
            exceed_t_pos = t_vals_perm>crit_t;
            exceed_t_pos = max(accumarray(nonzeros((cumsum(~exceed_t_pos)+1).*exceed_t_pos),1));
            if isempty(exceed_t_pos)
                exceed_t_pos = 0;
            end
            
            exceed_t_neg = t_vals_perm<-crit_t;
            exceed_t_neg = max(accumarray(nonzeros((cumsum(~exceed_t_neg)+1).*exceed_t_neg),1));
            if isempty(exceed_t_neg)
                exceed_t_neg = 0;
            end
            
            sig_t_maxClus(idx_perm) = max([exceed_t_pos exceed_t_neg]);
            
        end
        
        mm = beta(1,:);
        ss = sqrt(var_betas(1,:));
        
        mk_clusterStats(t_vals(1,:),sig_t_maxClus,true,true,19)
        
        figure, set(gcf,'color','white'),hold on
        if isempty(title_covariate)
            title('Covariate')
        else 
            title(title_covariate)
        end
        ylabel('Coefficient')
%         xlabel('lag (seconds)')
%         xticks([0:10:n_time])
%         xticklabels({(0:10:n_time)/100})  
        p = plot(mm,'Color',col{col_idx(idx_data)});
        if plot_se
            fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{col_idx(idx_data)},'EdgeAlpha',0,'FaceAlpha',0.2);
        end 
        plot(zeros(1,n_time),'--','Color',col{end-1});
        
        if diff_plot
%         ylabel('difference')                    
        if ~strcmp(which_data,'rest')
            xticks([1 6:5:n_time])
    %         xticklabels({(0:5:n_time)*100})
            xticklabels({(-5:5:n_time)*100})
            xlabel('time interval, 500ms  starting at (ms)')
        else
%                         xticks([1 11:10:n_time])
            xticks(linspace(1,n_time,9))
            xticklabels({(0:0.5:4)})
            xlabel('time interval, 0.5 minutes  starting at (minutes)')
        end
    else
%         ylabel('sequenceness')
        xlabel('lag (seconds)')
        xticks([0:10:n_time])
        xticklabels({(0:10:n_time)/100})
    end
                
    end


end