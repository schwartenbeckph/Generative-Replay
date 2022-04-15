% do RSA on timeseries, cf. Luyckx et al, eLife, 2019


function [RSM, EL_RSM, REL_RSM, Chunk_RSM, Graph_RSM, PIXEL_RSM, SIZE_RSM, STIM_shape, ... 
          RSM_sub1, RSM_sub2, PIXEL_RSM_sub1, SIZE_RSM_sub1, REL_RSM_sub1, PIXEL_RSM_sub2, SIZE_RSM_sub2, REL_RSM_sub2, beta_all] = do_RSA_timeseries_inference(scan_result_path,name_file,which_chan,lab_chunk,STIMS_all,plot_RSMs,which_part)

    if nargin<8
        plot_RSMs = false;
        graph = 0;
        which_part = [];
    end

    %% prelim %%
    map = [1 1 1; 
           linspace(0,1,10-1)' ones(10-1,1)*0.5    linspace(1,0,10-1)'];   
    map2 = [linspace(1,0,10)' linspace(1,0,10)'    linspace(1,0,10)']; % grey
    
    map3 = [linspace(0.8500,0,10)' linspace(0.50,0.4470,10)' linspace(0.0980,0.7410,10)'];
    
%     map3 = [linspace(0,1,10)' linspace(0,1,10)' linspace(1,0,10)'];
%     map3 = [linspace(0,0.8500,10)' linspace(0.4470,0.3250,10)' linspace(0.7410,0.0980,10)'];
%     map3 = [linspace(0,0.8500,10)' linspace(0.4470,0.50,10)' linspace(0.7410,0.0980,10)'];
%     map3 = [linspace(0,0.9290,10)' linspace(0.4470,0.6940,10)' linspace(0.7410,0.1250,10)'];

    load(fullfile(scan_result_path,name_file)) % this loads 'data' and 'stimlabel'
    
    if strcmp(which_part,'first_half')
        data = data(:,:,1:3*48);
        stimlabel = stimlabel(1:3*48);
        bricks_conn_trial = bricks_conn_trial(1:3*48,:);
        bricks_rel_trial = bricks_rel_trial(1:3*48,:);
        correct_trials_all = correct_trials_all(1:3*48);
        STIMS_all = STIMS_all(:,:,1:3*48);
    elseif strcmp(which_part,'second_half')
        data = data(:,:,3*48+1:end);
        stimlabel = stimlabel(3*48+1:end);
        bricks_conn_trial = bricks_conn_trial(3*48+1:end,:);
        bricks_rel_trial = bricks_rel_trial(3*48+1:end,:);
        correct_trials_all = correct_trials_all(3*48+1:end);
        STIMS_all = STIMS_all(:,:,3*48+1:end);        
    end
    
    %%%%% take specific channels? %%%%% 
    channel_names = char(channel_names');
    
    if strcmp(which_chan,'all')
        idx_chan = true(size(data,1),1);
    elseif strcmp(which_chan,'Occ')
        idx_chan = channel_names(:,3)=='O';
    elseif strcmp(which_chan,'Front')
        idx_chan = channel_names(:,3)=='F'; 
    elseif strcmp(which_chan,'Temp')
        idx_chan = channel_names(:,3)=='T'; 
    elseif strcmp(which_chan,'Par')
        idx_chan = channel_names(:,3)=='P'; 
    elseif strcmp(which_chan,'Cen')
        idx_chan = channel_names(:,3)=='C';         
    end
    
    data = data(idx_chan,:,:);    
    
    %%%%% baseline correct (optional) %%%%%
    data = baseline_correct_timeSeries(data,[41,50]); % 51 is is stim onset   
    
    %%%%% remove trials with broken recording %%%%%
    idx_Nonan     = squeeze(~isnan(data(1,51,:)));
    
    data = data(:,:,idx_Nonan);
    
    labStm    = zeros(1,length(stimlabel));
    allconds  = unique(stimlabel);
    num_label = 1:length(allconds);

    for idx_num_label=1:max(num_label)
        labStm(strcmp(stimlabel,allconds{idx_num_label})) = idx_num_label;
    end
    labStm        = labStm(idx_Nonan);
    
    bricks_conn_trial = bricks_conn_trial(idx_Nonan,:);
    bricks_rel_trial  = bricks_rel_trial(idx_Nonan,:);
    
    [n_chan,n_time,n_trial] = size(data);
    
    %%%%% z-score the data across trials (cf., Luyckx et al, eLife, 2019) %%%%%
%     for idx_chan=1:n_chan
%         for idx_time=1:n_time
%             data(idx_chan,idx_time,:) = zscore(data(idx_chan,idx_time,:));
%         end
%     end
    
    %%%%% make design matrix %%%%%
    X = zeros(n_trial,max(labStm));
    
    for idx_stim=1:max(labStm)
        X(:,idx_stim) = labStm==idx_stim;
    end
    
    %% make theoretical RSMs %%
    
    STIMS_all = STIMS_all(:,:,idx_Nonan);
    
    El_stim    = zeros(max(labStm),3);
    ElR_stim   = zeros(max(labStm),4);
    STIM_shape = zeros(6,6,max(labStm));
    for idx_stim=1:max(labStm)
        El_stim(idx_stim,:) = sort(unique(bricks_conn_trial(labStm==idx_stim,:),'rows'));
        ElR_stim(idx_stim,:) = unique(bricks_rel_trial(labStm==idx_stim,:),'rows');
        
        STIMS_idx = STIMS_all(:,:,labStm==idx_stim);
%         mk_STIM_plot(STIMS_idx)
        
        STIM_shape(:,:,idx_stim)  = STIMS_idx(:,:,1);
        
        STIM_shape_cell{idx_stim} = STIMS_idx(:,:,1);
    end
    
%     mk_STIM_plot(STIM_shape) 

    % 0) visual similarities (pixel and size)
    [diff_heightwidth, Vis_similarity_move] = mk_vis_sim(1:max(labStm),STIM_shape_cell);    
    
    PIXEL_RSM = Vis_similarity_move;
    SIZE_RSM  = 1-diff_heightwidth;
    
    if plot_RSMs
        mk_RSA_plot([1:12],PIXEL_RSM,STIM_shape_cell)   
        colormap(map2)
%         colormap(map3)
        mk_RSA_plot([1:12],SIZE_RSM,STIM_shape_cell)   
        colormap(map2)
%         colormap(map3)
    end
    
    PIXEL_RSM_mat = PIXEL_RSM;
    SIZE_RSM_mat  = SIZE_RSM;
    
    PIXEL_RSM = PIXEL_RSM(tril(true(length(PIXEL_RSM)),-1));
    SIZE_RSM  = SIZE_RSM(tril(true(length(SIZE_RSM)),-1));
    
    % i) same elements
    [~,~,idx_el] = unique(El_stim,'rows');
    
    EL_RSM = zeros(max(labStm));
    for idx_row=1:max(labStm)
        for idx_col=1:max(labStm)
            EL_RSM(idx_row,idx_col) = idx_el(idx_row) == idx_el(idx_col);
        end
    end
    
%     EL_RSM = sqrt(idx_el*idx_el')==floor(sqrt(idx_el*idx_el'));
    
    if plot_RSMs
        mk_RSA_plot([1:12],EL_RSM,STIM_shape_cell)   
        colormap(map3)
    end
    
    EL_RSM = EL_RSM(tril(true(length(EL_RSM)),-1));
    
    % ii) element in relation overlap
    REL_RSM = mk_ELinREL(1:max(labStm),ElR_stim);
    
    if plot_RSMs
        mk_RSA_plot([1:12],REL_RSM,STIM_shape_cell)
        colormap(map2)
%         colormap(map3)
    end
    
    REL_RSM_mat = REL_RSM;
    REL_RSM     = REL_RSM(tril(true(length(REL_RSM)),-1));
    
    % iii) chunk overlap
    Chunk_stim = zeros(max(labStm),2);
    for idx_stim=1:max(labStm)
        Chunk_stim(idx_stim,:) = [find(ismember(lab_chunk(:,1:2),El_stim(idx_stim,1:2),'rows')); ...
                                  find(ismember(lab_chunk(:,1:2),El_stim(idx_stim,2:3),'rows')); ...
                                  find(ismember(lab_chunk(:,1:2),fliplr(El_stim(idx_stim,1:2)),'rows')); ...
                                  find(ismember(lab_chunk(:,1:2),fliplr(El_stim(idx_stim,2:3)),'rows'))]';
    end
    
    Chunk_RSM = zeros(max(labStm),max(labStm));
    for idx_col=1:max(labStm)
       for idx_row=1:max(labStm)
           Chunk_RSM(idx_row,idx_col) = mean(Chunk_stim(idx_row,:) == Chunk_stim(idx_col,:));
       end
    end
    
    if plot_RSMs
        mk_RSA_plot([1:12],Chunk_RSM,STIM_shape_cell)
        colormap(map2)
    end
    
    Chunk_RSM = Chunk_RSM(tril(true(length(Chunk_RSM)),-1));
    
    % iv) graph overlap
    labStm_struct = nan(max(labStm),1);
    % obtain 'graph structures'
    labStm_struct(ElR_stim(:,1)==ElR_stim(:,2)) = 1;
    labStm_struct(ElR_stim(:,1)==ElR_stim(:,4)) = 2;
    labStm_struct(ElR_stim(:,2)==ElR_stim(:,3)) = 3;
    labStm_struct(ElR_stim(:,3)==ElR_stim(:,4)) = 4;
    
    Graph_RSM = zeros(max(labStm));
    for idx_row=1:max(labStm)
        for idx_col=1:max(labStm)
            Graph_RSM(idx_row,idx_col) = labStm_struct(idx_row) == labStm_struct(idx_col);
        end
    end
    
%     Graph_RSM = sqrt(labStm_struct*labStm_struct')==floor(sqrt(labStm_struct*labStm_struct'));
    
    if plot_RSMs
        mk_RSA_plot([1:12],Graph_RSM,STIM_shape_cell)
        colormap(map2)
    end
    
    Graph_RSM = Graph_RSM(tril(true(length(Graph_RSM)),-1));
    
    %% find uncorrelatated elements in relations and pixel/size
    % quick and dirty:
    diff_rel_pixel = REL_RSM_mat-PIXEL_RSM_mat;
    diff_rel_size  = REL_RSM_mat-SIZE_RSM_mat;
    
    [~,idx_pixel] = sort(sum(abs(diff_rel_pixel),2));
    [~,idx_size]  = sort(sum(abs(diff_rel_size),2));

    idx_sub1 = idx_pixel(7:end)'; % part of data with small correlation
    idx_sub2 = idx_pixel(1:6)'; % part of data with large correlation
    
    [diff_heightwidth_sub1, Vis_similarity_move_sub1] = mk_vis_sim(idx_sub1,STIM_shape_cell);        
    PIXEL_RSM_sub1 = Vis_similarity_move_sub1;
    SIZE_RSM_sub1  = 1-diff_heightwidth_sub1;
    
    % exact:
%     if ~exist(fullfile('D:\StimConstrMEG\Results_Scanning\Summary_Data',['idx_graph_',num2str(graph),'.mat']),'file')
%         all_idx = nchoosek(1:12,6);
%         current_min = 10;
%         idx_sub_min = zeros(1,6);
% 
%         for idx_nchoosek=1:size(all_idx,1)
% 
%             % choose six stims
%             idx_sub1 = all_idx(idx_nchoosek,:);        
% 
%             % obtain visual info
%             [diff_heightwidth_sub1, Vis_similarity_move_sub1] = mk_vis_sim(idx_sub1,STIM_shape_cell);        
%             PIXEL_RSM_sub1 = Vis_similarity_move_sub1;
%             SIZE_RSM_sub1  = 1-diff_heightwidth_sub1;
% 
%             % obtain stats:
%             if current_min>mk_meanCorr_PS(PIXEL_RSM_sub1(tril(true(length(PIXEL_RSM_sub1)),-1)))
%                 idx_sub_min = idx_sub1;
%                 current_min = mk_meanCorr_PS(PIXEL_RSM_sub1(tril(true(length(PIXEL_RSM_sub1)),-1)));
%             end
% %             if current_min>(mk_meanCorr_PS(PIXEL_RSM_sub1(tril(true(length(PIXEL_RSM_sub1)),-1)))+SIZE_RSM_sub1(tril(true(length(SIZE_RSM_sub1)),-1)))
% %                 idx_sub_min = idx_sub1;
% %                 current_min = mk_meanCorr_PS(PIXEL_RSM_sub1(tril(true(length(PIXEL_RSM_sub1)),-1)))+SIZE_RSM_sub1(tril(true(length(SIZE_RSM_sub1)),-1));
% %             end
% 
%             fprintf('Idx %d of %d done.\n',idx_nchoosek,size(all_idx,1))
% 
%         end
% 
%         idx_sub1 = idx_sub_min;
%         idx_sub2 = setdiff(1:12,idx_sub_min);
%         
%         save(fullfile('D:\StimConstrMEG\Results_Scanning\Summary_Data',['idx_graph_',num2str(graph),'.mat']),...
%              'idx_sub1','idx_sub2','-v7.3')
%         
%     else
%         load(fullfile('D:\StimConstrMEG\Results_Scanning\Summary_Data',['idx_graph_',num2str(graph),'mat']))
%     end
    
    REL_RSM_sub1 = mk_ELinREL(idx_sub1,ElR_stim);
    
    if plot_RSMs
        mk_RSA_plot(idx_sub1,PIXEL_RSM_sub1,STIM_shape_cell)   
        colormap(map2)
        mk_RSA_plot(idx_sub1,SIZE_RSM_sub1,STIM_shape_cell)   
        colormap(map2)
        mk_RSA_plot(idx_sub1,REL_RSM_sub1,STIM_shape_cell)
        colormap(map2)
    end
    
    PIXEL_RSM_sub1 = PIXEL_RSM_sub1(tril(true(length(PIXEL_RSM_sub1)),-1));
    SIZE_RSM_sub1  = SIZE_RSM_sub1(tril(true(length(SIZE_RSM_sub1)),-1));
    REL_RSM_sub1   = REL_RSM_sub1(tril(true(length(REL_RSM_sub1)),-1));
    
    [diff_heightwidth_sub2, Vis_similarity_move_sub2] = mk_vis_sim(idx_sub2,STIM_shape_cell);        
    PIXEL_RSM_sub2 = Vis_similarity_move_sub2;
    SIZE_RSM_sub2  = 1-diff_heightwidth_sub2;
    
    REL_RSM_sub2 = mk_ELinREL(idx_sub2,ElR_stim);
    
    if plot_RSMs
        mk_RSA_plot(idx_sub2,PIXEL_RSM_sub2,STIM_shape_cell)   
        colormap(map2)
        mk_RSA_plot(idx_sub2,SIZE_RSM_sub2,STIM_shape_cell)   
        colormap(map2)
        mk_RSA_plot(idx_sub2,REL_RSM_sub2,STIM_shape_cell)
        colormap(map2)
    end    
    
    PIXEL_RSM_sub2 = PIXEL_RSM_sub2(tril(true(length(PIXEL_RSM_sub2)),-1));
    SIZE_RSM_sub2  = SIZE_RSM_sub2(tril(true(length(SIZE_RSM_sub2)),-1));
    REL_RSM_sub2   = REL_RSM_sub2(tril(true(length(REL_RSM_sub2)),-1));
    
    all_reg  = [PIXEL_RSM,SIZE_RSM,REL_RSM];
    reg_sub1 = [PIXEL_RSM_sub1,SIZE_RSM_sub1,REL_RSM_sub1];
    reg_sub2 = [PIXEL_RSM_sub2,SIZE_RSM_sub2,REL_RSM_sub2];
    
    r      = corr(all_reg);
    r_sub1 = corr(reg_sub1);
    r_sub2 = corr(reg_sub2);
    
    fprintf('Overall corr is %f %f %f.\n',r(2,1),r(3,1),r(3,2))
    fprintf('Reduced corr is %f %f %f.\n',r_sub1(2,1),r_sub1(3,1),r_sub1(3,2))
    fprintf('Larger corr is %f %f %f.\n',r_sub2(2,1),r_sub2(3,1),r_sub2(3,2))

    %% make empirical RSMs %%

    [RSM, beta_all] = mk_RSA_timeseries(labStm,n_time,data,X,idx_chan_use,false);
    
    RSM_sub1 = mk_RSA_timeseries(labStm(ismember(labStm,idx_sub1)),n_time,data(:,:,ismember(labStm,idx_sub1)),X(ismember(labStm,idx_sub1),idx_sub1),idx_chan_use,false);
    
    RSM_sub2 = mk_RSA_timeseries(labStm(ismember(labStm,idx_sub2)),n_time,data(:,:,ismember(labStm,idx_sub2)),X(ismember(labStm,idx_sub2),idx_sub2),idx_chan_use,false);

end