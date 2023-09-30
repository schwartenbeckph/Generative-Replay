
function [] = plot_timeseries_PS(data,Y_name,vert_idx,legend_text,x_Tick,title_lab,do_perm,one_sided,perc_data,col_idx)

    col   = {[0, 0.4470, 0.7410], ...       % blue
             [0.4940, 0.1840, 0.5560], ...  % purple
             [0.4660, 0.6740, 0.1880], ...  % green
             [0.9350, 0.1780, 0.2840], ...  % red
             [0.3010, 0.7450, 0.9330], ...  % cyan
             [0.9290, 0.6940, 0.1250], ...  % yellow
             [0, 0, 0], ...                 % black
             [0.6, 0.6, 0.6]};              % grey 

    if nargin<10
        col_idx = 1:6;
    end
    
    if nargin<9
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<8
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<7
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<6
        title_lab = [];
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end
    
    if nargin<5
        x_Tick = [];
        title_lab = [];
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<4
        legend_text = [];
        x_Tick = [];
        title_lab = [];
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<3
        vert_idx = [];
        legend_text = [];
        x_Tick = [];
        title_lab = [];
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if nargin<2
        vert_idx = [];
        legend_text = [];
        x_Tick = [];
        Y_name = 'measure';
        title_lab = [];
        do_perm = false;
        one_sided = false;
        perc_data = [97.5 2.5];
        col_idx = 1:6;
    end

    if ~do_perm
        [n_sub,n_time,n_data] = size(data);
    else
        [n_sub,n_time,n_perm,n_data] = size(data);
    end
    
    col_idx = col_idx(1:n_data);
    
    if do_perm
        for idx_dat=1:n_data
            [up_sess(:,idx_dat), low_sess(:,idx_dat), emp_sess(:,idx_dat)] = sign_flip_AcrossTime_PS([],[],data(:,:,:,idx_dat),perc_data);
        end
    end
    
    if exist('n_perm','var')
       data = squeeze(data(:,:,1,:)); 
    end    
         
    y_max = max(max(max(squeeze(nanmean(data,1)))));
    if exist('up_sess','var')
        y_max = max(y_max,max(max(up_sess)));
    end
    y_max = round(y_max*1.5,2);
%     y_max = round(y_max*2,2);
%     y_max = round(y_max*3,2);

    y_min = min(min(min(squeeze(nanmean(data,1)))));
    if exist('low_sess','var')
        y_min = min(y_min,min(min(low_sess)));
    end
    y_min = round(y_min*1.5,2);
%     y_min = round(y_min*2,2);
%     y_min = round(y_min*3,2);    

    figure, set(gcf,'color','white'), hold on
    
    for idx_dat=1:n_data
        
        mm = nanmean(data(:,:,idx_dat),1);
        ss = nanstd(data(:,:,idx_dat),0,1)./sqrt(n_sub);
        
        p(idx_dat) = plot(mm,'Color',col{col_idx(idx_dat)},'LineWidth',3);
        fill([1:n_time, fliplr(1:n_time)],[mm+ss, fliplr(mm-ss)], col{col_idx(idx_dat)},'EdgeAlpha',0,'FaceAlpha',0.2);
        
        if n_data==1
            fill([1:n_time, fliplr(1:n_time)],[up_sess(:,idx_dat)', fliplr(low_sess(:,idx_dat))'], col{8},'EdgeAlpha',0,'FaceAlpha',0.2);
        end
        if one_sided
            if n_data==1
                plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(:,idx_dat),95),'--k','MarkerSize',6) % one sided test
            else
                plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(:,idx_dat),95),'--','Color',col{col_idx(idx_dat)},'MarkerSize',6) % one sided test
            end
        else
            if n_data==1
                plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(:,idx_dat),97.5),'--k','MarkerSize',6) % one sided test
                plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess(:,idx_dat),2.5),'--k','MarkerSize',6) % one sided test
            else
                plot(1:n_time,ones(length(1:n_time),1)*prctile(up_sess(:,idx_dat),97.5),'--','Color',col{col_idx(idx_dat)},'MarkerSize',6) % one sided test
                plot(1:n_time,ones(length(1:n_time),1)*prctile(low_sess(:,idx_dat),2.5),'--','Color',col{col_idx(idx_dat)},'MarkerSize',6) % one sided test
            end
        end
        
    end
    
    ylim([min(y_min),max(y_max)])
    
    if ~isempty(vert_idx)
        for idx_vert_idx = vert_idx
            line([idx_vert_idx idx_vert_idx],get(gca, 'ylim'),'Color','k');
        end
    end
    
    line(get(gca, 'xlim'),[0 0],'Color','k','LineStyle','--');

    xlabel('time')
    xticks(0:50:n_time)
    if ~isempty(x_Tick)
        x_Labels = linspace(x_Tick(1),x_Tick(2),length(0:50:n_time));
        xticklabels({x_Labels})
    end
    ylabel(Y_name)
    
    if ~isempty(legend_text)
        legend([p],legend_text)
    end
    
    if ~isempty(title_lab)
        title(title_lab)
    end
    
    hold off

end