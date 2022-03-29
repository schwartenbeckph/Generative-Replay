
function [sig_where] = mk_clusterStats(input_data,sig_t_maxClus,do_plot,data_is_t,df,p_crit,diff_data,one_sided,cols_use,title_input,legend_input,y_axis_Input)

    sig_where{1} = [];
    sig_where{2} = [];

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
         
    [n_sub, n_time, n_data] = size(input_data);

     if nargin<4
        data_is_t = false;
        df = size(input_data,1)-1; % implies all your data must have same number of subjects
        p_crit = 0.05;
        diff_data = false;
        one_sided = false;
        cols_use = 1:n_data;
        title_input = [];
        legend_input = [];
        y_axis_Input = [];
     elseif nargin<9
        cols_use = 1:n_data;
        title_input = [];
        legend_input = [];
        y_axis_Input = [];
    elseif nargin<12
        y_axis_Input = [];
    end
    
    if do_plot
        figure, set(gcf,'color','white')
        hold on
        
        y_ax = nanmean(input_data,1);
        y_ax = max([max(max(y_ax)),max(max(-y_ax))]);
    end
    
    for idx_data=1:n_data

        % now find entries that exceed cluster thresh
        cluster_bound = prctile(sig_t_maxClus(1,2:end,idx_data),95);
        cluster_bound = round(cluster_bound);
        t_diff_thresh = zeros(1,size(input_data(:,:,idx_data),2)+2*cluster_bound);
        
%         fprintf('cluster bound is %d.\n',cluster_bound)

        if data_is_t
            t_diff = input_data(:,:,idx_data);
        else
            t_diff = mean(input_data(:,:,idx_data))./(std(input_data(:,:,idx_data))/sqrt(size(input_data(:,:,idx_data),1)));
        end

        if strcmp(one_sided,'one_sided_neg')

            t_diff_thresh_pos = false(size(t_diff));
            t_diff_thresh_neg = t_diff<=tinv(p_crit,df);

        elseif strcmp(one_sided,'one_sided_pos')

            t_diff_thresh_pos = t_diff>=tinv(1-p_crit,df);
            t_diff_thresh_neg = false(size(t_diff));

        else

            t_diff_thresh_pos = t_diff>=tinv(1-p_crit/2,df);
            t_diff_thresh_neg = t_diff<=tinv(p_crit/2,df);

        end

        t_diff_thresh_pos = [zeros(1,cluster_bound) t_diff_thresh_pos zeros(1,cluster_bound)];

        t_diff_thresh_neg = [zeros(1,cluster_bound) t_diff_thresh_neg zeros(1,cluster_bound)];

        for idx_nonzero=find(t_diff_thresh_pos)
            data_temp = t_diff_thresh_pos(idx_nonzero-cluster_bound:idx_nonzero+cluster_bound);
            t_diff_thresh(idx_nonzero) = max(accumarray(nonzeros((cumsum(~data_temp)+1).*data_temp),1));
        end

        for idx_nonzero=find(t_diff_thresh_neg)
            data_temp = t_diff_thresh_neg(idx_nonzero-cluster_bound:idx_nonzero+cluster_bound);
            t_diff_thresh(idx_nonzero) = -max(accumarray(nonzeros((cumsum(~data_temp)+1).*data_temp),1));
        end

        t_diff_thresh = t_diff_thresh(cluster_bound+1:end-cluster_bound);       

        if do_plot

            mm = nanmean(input_data(:,:,idx_data),1);
            ss = nanstd(input_data(:,:,idx_data),0,1)./sqrt(n_sub);
            p(idx_data) = plot(mm,'Color',col{cols_use(idx_data)},'LineWidth',3);
            fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{cols_use(idx_data)},'EdgeAlpha',0,'FaceAlpha',0.2);
            plot(zeros(1,length(t_diff)),'--k')


            if ~isempty(find(t_diff_thresh>=cluster_bound))
                if n_data==1
                    plot(find(t_diff_thresh>=cluster_bound),(y_ax)*ones(size(find(t_diff_thresh>=cluster_bound))),'.','Color',col{cols_use(idx_data)},'LineWidth',5)
                    
                    sig_where{1} = find(t_diff_thresh>=cluster_bound);
                else
                    plot(find(t_diff_thresh>=cluster_bound),(y_ax*(1+idx_data/20))*ones(size(find(t_diff_thresh>=cluster_bound))),'.','Color',col{cols_use(idx_data)},'LineWidth',5)
                    
                    sig_where{1}{idx_data} = find(t_diff_thresh>=cluster_bound);
                end
            end
            if ~isempty(find(t_diff_thresh<=-cluster_bound))
                if n_data==1
                    plot(find(t_diff_thresh<=-cluster_bound),(-y_ax)*ones(size(find(t_diff_thresh<=-cluster_bound))),'.','Color',col{cols_use(idx_data)},'LineWidth',5)
                    
                    sig_where{2} = find(t_diff_thresh<=-cluster_bound);
                else
                    plot(find(t_diff_thresh<=-cluster_bound),(-y_ax*(idx_data/10))*ones(size(find(t_diff_thresh<=-cluster_bound))),'.','Color', col{cols_use(idx_data)},'LineWidth',5)
                    
                    sig_where{2}{idx_data} = find(t_diff_thresh>=cluster_bound);
                end
            end

        end
    
    end
    
    if do_plot
        
        if ~isempty(y_axis_Input)
            ylim([y_axis_Input(1),y_axis_Input(2)])
        else
            ylim([-y_ax*2,y_ax*2])
        end
        
        if diff_data
            ylabel('Sequenceness')
            xlabel('time interval, 1000ms  starting at (ms)')
            xticks([51 151 251])
            xticklabels({'[0,1000]' '[1000,2000]' '[2000,3000]'})
            if isempty(title_input)
                title('Difference Seq from to-be-inferred and stable building blocks')
            else
                title(title_input)
            end
            if ~isempty(legend_input)
                legend(p,legend_input)
            end
        else
            ylabel('sequenceness')
            xlabel('lag (seconds)')
            xticks([0:10:size(input_data,2)])
            xticklabels({(0:10:size(input_data,2))/100})
        end
        
        hold off
    end


end