% obtain brick info of subject based on observed stims
function [bricks_conn,bricks_rel,bricks_conn_trial,bricks_rel_trial,bricks_conn_Q,bricks_rel_Q,correct_trials_all,rt_all,stims_all,STIMS_all,relation_catch_all] = mk_bricks_sub(behav_folder)

    behav_files    = dir(fullfile(behav_folder,'T*.mat'));

    correct_trials_all = [];

    STIMS_all = [];
    rt_all    = [];
    
    stims_all = [];
    
    stim_catch_all     = [];
    relation_catch_all = [];

    for idx_sess = 1:length(behav_files)

        load(fullfile(behav_folder,behav_files(idx_sess).name))

        correct_trials_all = [correct_trials_all; res.behav.correct'];

        STIMS_all = cat(3,STIMS_all,res.behav.SOLUTIONS_BUILT);

        rt_all = [rt_all; res.behav.rt'];
        
        stim_catch_all = [stim_catch_all; res.behav.stim_catch];
        
        relation_catch_all = [relation_catch_all; res.behav.relation_catch'];
        
        stims_all = [stims_all; res.order];

    end  

    %%%%% now obtain unique stims and then get brick ordering %%%%%
    [n,m,p] = size(STIMS_all);
    a = reshape(STIMS_all,n,[],1);
    b = reshape(a(:),n*m,[])';
    c = unique(b,'rows','stable')';
    STIMS_unique = reshape(c,n,m,[]);

    % obtain summary for all unique stims shown
    
    % there are two ways of ordering the stims:
    bricks_conn = zeros(size(STIMS_unique,3),3); % Element connected to middle via besideness | middle Element | Element connected to middle via ontopness
    bricks_rel  = zeros(size(STIMS_unique,3),4); % left element | ontop element | right element | below element
    
    for idx_stim=1:size(STIMS_unique,3)

        STIMS_obj = STIMS_unique(:,:,idx_stim);

        bricks = unique(STIMS_obj); bricks(1) = []; % get rid of zero

        PART_no1 = STIMS_obj(:,:,1); PART_no1(PART_no1==bricks(1)) = 0;

        PART_no2 = STIMS_obj(:,:,1); PART_no2(PART_no2==bricks(2)) = 0;

        PART_no3 = STIMS_obj(:,:,1); PART_no3(PART_no3==bricks(3)) = 0;

        % order bricks such that first two are connected via besideness and second two via ontopess, and middle column is the stimulus in middle
        % isolate parts that contain brick 1 | brick 2 | brick 3, see if they contain ontopness; besideness
        bricks_order = [mk_ontopness(PART_no3)+mk_ontopness(PART_no2)  mk_ontopness(PART_no1)+mk_ontopness(PART_no3)   mk_ontopness(PART_no1)+mk_ontopness(PART_no2);
                        mk_besideness(PART_no3)+mk_besideness(PART_no2) mk_besideness(PART_no1)+mk_besideness(PART_no3) mk_besideness(PART_no1)+mk_besideness(PART_no2)];

        bricks_order = [find(~bricks_order(1,:)&bricks_order(2,:)) find(bricks_order(1,:)&bricks_order(2,:)) find(bricks_order(1,:)&~bricks_order(2,:))];

        bricks_conn(idx_stim,:) = bricks(bricks_order)';  

        if mk_ontopness(PART_no1)
            [~,~,bricks_rel(idx_stim,2),bricks_rel(idx_stim,4)] = mk_ontopness(PART_no1);
        elseif mk_ontopness(PART_no2)
            [~,~,bricks_rel(idx_stim,2),bricks_rel(idx_stim,4)] = mk_ontopness(PART_no2);
        elseif mk_ontopness(PART_no3)
            [~,~,bricks_rel(idx_stim,2),bricks_rel(idx_stim,4)] = mk_ontopness(PART_no3);
        end

        if mk_besideness(PART_no1)
            [~,~,bricks_rel(idx_stim,1),bricks_rel(idx_stim,3)] = mk_besideness(PART_no1);
        elseif mk_besideness(PART_no2)
            [~,~,bricks_rel(idx_stim,1),bricks_rel(idx_stim,3)] = mk_besideness(PART_no2);
        elseif mk_besideness(PART_no3)
            [~,~,bricks_rel(idx_stim,1),bricks_rel(idx_stim,3)] = mk_besideness(PART_no3);
        end            

    end
    
    % obtain summary for all stims per trial shown
    
    % there are two ways of ordering the stims:
    bricks_conn_trial = zeros(size(STIMS_all,3),3); % Element connected to middle via besideness | middle Element | Element connected to middle via ontopness
    bricks_rel_trial  = zeros(size(STIMS_all,3),4); % left element | ontop element | right element | below element
    
    bricks_conn_Q = zeros(size(STIMS_all,3),3); % Element connected to middle via besideness | middle Element | Element connected to middle via ontopness
    bricks_rel_Q  = zeros(size(STIMS_all,3),4); % left element | ontop element | right element | below element
    
    for idx_stim=1:size(STIMS_all,3)

        STIMS_obj = STIMS_all(:,:,idx_stim);

        bricks = unique(STIMS_obj); bricks(1) = []; % get rid of zero

        PART_no1 = STIMS_obj(:,:,1); PART_no1(PART_no1==bricks(1)) = 0;

        PART_no2 = STIMS_obj(:,:,1); PART_no2(PART_no2==bricks(2)) = 0;

        PART_no3 = STIMS_obj(:,:,1); PART_no3(PART_no3==bricks(3)) = 0;

        % order bricks such that first two are connected via besideness and second two via ontopess, and middle column is the stimulus in middle
        % isolate parts that contain brick 1 | brick 2 | brick 3, see if they contain ontopness; besideness
        bricks_order = [mk_ontopness(PART_no3)+mk_ontopness(PART_no2)  mk_ontopness(PART_no1)+mk_ontopness(PART_no3)   mk_ontopness(PART_no1)+mk_ontopness(PART_no2);
                        mk_besideness(PART_no3)+mk_besideness(PART_no2) mk_besideness(PART_no1)+mk_besideness(PART_no3) mk_besideness(PART_no1)+mk_besideness(PART_no2)];

        bricks_order = [find(~bricks_order(1,:)&bricks_order(2,:)) find(bricks_order(1,:)&bricks_order(2,:)) find(bricks_order(1,:)&~bricks_order(2,:))];

        bricks_conn_trial(idx_stim,:) = bricks(bricks_order)';  

        if mk_ontopness(PART_no1)
            [~,~,bricks_rel_trial(idx_stim,2),bricks_rel_trial(idx_stim,4)] = mk_ontopness(PART_no1);
        elseif mk_ontopness(PART_no2)
            [~,~,bricks_rel_trial(idx_stim,2),bricks_rel_trial(idx_stim,4)] = mk_ontopness(PART_no2);
        elseif mk_ontopness(PART_no3)
            [~,~,bricks_rel_trial(idx_stim,2),bricks_rel_trial(idx_stim,4)] = mk_ontopness(PART_no3);
        end

        if mk_besideness(PART_no1)
            [~,~,bricks_rel_trial(idx_stim,1),bricks_rel_trial(idx_stim,3)] = mk_besideness(PART_no1);
        elseif mk_besideness(PART_no2)
            [~,~,bricks_rel_trial(idx_stim,1),bricks_rel_trial(idx_stim,3)] = mk_besideness(PART_no2);
        elseif mk_besideness(PART_no3)
            [~,~,bricks_rel_trial(idx_stim,1),bricks_rel_trial(idx_stim,3)] = mk_besideness(PART_no3);
        end 
        
        bricks_conn_Q(idx_stim,:) = bricks_conn_trial(idx_stim,:).*ismember(bricks_conn_trial(idx_stim,:),stim_catch_all(idx_stim,:));
        bricks_rel_Q(idx_stim,:)  = bricks_rel_trial(idx_stim,:).*ismember(bricks_rel_trial(idx_stim,:),stim_catch_all(idx_stim,:));
        
        % make sure bricks_rel_Q only contains question relevant relation positions
        if relation_catch_all(idx_stim)==1 || relation_catch_all(idx_stim)==3
            bricks_rel_Q(idx_stim,2) = 0;
            bricks_rel_Q(idx_stim,4) = 0;
        elseif relation_catch_all(idx_stim)==2 || relation_catch_all(idx_stim)==4
            bricks_rel_Q(idx_stim,1) = 0;
            bricks_rel_Q(idx_stim,3) = 0;
        end
         
    end    

end