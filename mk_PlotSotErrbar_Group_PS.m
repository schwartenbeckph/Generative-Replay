

function [] = mk_PlotSotErrbar_Group_PS(input_data,horiz_line,input_ylabel,input_xlabel,input_title,input_XTick,input_h,col_idx,input_legend,input_y_lim_ax,input_x_lim_ax,plot_chance)

%     col   = {[0, 0.4470, 0.7410], ...       % blue
%              [0.4940, 0.1840, 0.5560], ...  % purple
%              [0.4660, 0.6740, 0.1880], ...  % green
%              [0.9350, 0.1780, 0.2840], ...  % red
%              [0.3010, 0.7450, 0.9330], ...  % cyan
%              [0.9290, 0.6940, 0.1250], ...  % yellow
%              [0, 0, 0], ...                 % black
%              [0.6, 0.6, 0.6]};              % grey 

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


    if nargin==1
        horiz_line = false;
        input_ylabel = '';
        input_xlabel = '';
        input_title  = '';
        input_XTick  = [];
        input_h      = [];
        col_idx      = 1;
        input_legend = [];
        input_y_lim_ax = [];
        input_x_lim_ax = [];
        plot_chance = [];
    end
    
    figure
    set(gcf,'color','white')
    hold on
    
    if size(input_data,3)==1
        x_Ticks = 1:size(input_data,2);
    else
        x_Ticks = repmat(1:size(input_data,2),size(input_data,3),1) + linspace(-0.05*size(input_data,3),0.05*size(input_data,3),size(input_data,3))';
    end
    
    for idx_dat=1:size(input_data,3)
    
        p(idx_dat) = plot(x_Ticks(idx_dat,:),mean(input_data(:,:,idx_dat)),'.','Color',col{col_idx(idx_dat)},'MarkerSize',20);

        er = errorbar(x_Ticks(idx_dat,:),mean(input_data(:,:,idx_dat)),std(input_data(:,:,idx_dat),0,1)/sqrt(size(input_data(:,:,idx_dat),1)),'Color',col{col_idx(idx_dat)},'Linewidth',3);    
        er.LineStyle = 'none';
    
    end
    
    if isempty(input_y_lim_ax)
        y_lim_ax = [floor(min(min(mean(input_data))))-0.33*min(min(mean(input_data))), ...
                    ceil(max(max(mean(input_data))))+0.33*max(max(mean(input_data)))];           
    else
        y_lim_ax = input_y_lim_ax;
    end
    if isempty(input_x_lim_ax)
        x_lim_ax = [0,size(input_data,2)+1];
    else
        x_lim_ax = input_x_lim_ax;
    end
    xlim(x_lim_ax)
    ylim(y_lim_ax)

    if horiz_line
        line(get(gca, 'xlim'),[0 0],'Color','k','LineStyle','--');
    end

    ylabel(input_ylabel)
    xlabel(input_xlabel)
    title(input_title)

    if strcmp(input_XTick,'NoX')
        set(gca, 'XTick', '')
    elseif ~isempty(input_XTick)
        set(gca, 'XTick', [1:size(input_data,2)])
        set(gca, 'XTickLabel', input_XTick)
    else
        set(gca, 'XTick', [1:size(input_data,2)])
    end

    if ~isempty(input_h)
%         y_sig = ones(1,sum(input_h))*ceil(max(max(mean(input_data))));
        if size(input_h,3)==1
            y_sig = ones(1,sum(input_h))*y_lim_ax(2)*.95;

            x_sig = 1:size(input_data,2);
            x_sig = x_sig(logical(input_h));
            
            plot(x_sig,y_sig,'*k','MarkerSize',10,'Linewidth',2)
%             plot(x_sig,y_sig,'*','Color',col{col_idx})
        else
            y_sig = ones(1,sum(input_h(:)))*y_lim_ax(2)*.95;
            
            h = squeeze(input_h)';
            
            x_sig = x_Ticks(:);
            x_sig = x_sig(logical(h));
            
            col_sig = repmat(col_idx,size(input_h,2),1)';
            col_sig = col_sig(logical(h));
%             x_sig = [x_Ticks(1,:) x_Ticks(2,:)];
%             x_sig = [x_Ticks(:,1) x_Ticks(:,2)];
%             x_sig = x_sig(logical(input_h(:)));

            for idx_sig=1:length(x_sig)
                plot(x_sig(idx_sig),y_sig(idx_sig),'*','Color',col{col_sig(idx_sig)},'MarkerSize',10,'Linewidth',2)
            end
        end

%         plot(x_sig,y_sig,'*k')
%         plot(x_sig,y_sig,'*','Color',col{col_sig})
    end
    
    if ~isempty(input_legend)
        legend(p,input_legend)
    end
    
    if ~isempty(plot_chance)
        plot(x_lim_ax(1):x_lim_ax(2),ones(1,x_lim_ax(2)+1)*plot_chance,'--k')
    end
    
    hold off
    
end