% Function to make visual similarity RSM

function [diff_heightwidth, Vis_similarity_move, diff_height, diff_width] = mk_vis_sim(stim_Indices,STIM)
    
    diff_heightwidth  = nan(length(stim_Indices));
    Vis_similarity_move = nan(length(stim_Indices));

    for idx_one=1:length(stim_Indices)

        for idx_two=1:length(stim_Indices)

            comp_one_size = size(mk_crop(STIM{stim_Indices(idx_one)}));
            comp_two_size = size(mk_crop(STIM{stim_Indices(idx_two)}));

            diff_height = dist(comp_one_size(1),comp_two_size(1));
            diff_width  = dist(comp_one_size(2),comp_two_size(2));

            diff_heightwidth(idx_one,idx_two)  = dist(comp_one_size(1),comp_two_size(1)) + dist(comp_one_size(2),comp_two_size(2));

            comp_one = find(mk_recentre(STIM{stim_Indices(idx_one)}))';
            comp_two = find(mk_recentre(STIM{stim_Indices(idx_two)}))';
            Vis_similarity_move(idx_one,idx_two) = mk_visual_similarity_move(comp_one,comp_two);

%             fprintf('#### Index one %d, index two %d done ####\n',idx_one,idx_two)

        end

    end     

    diff_heightwidth = diff_heightwidth/max(max(diff_heightwidth));  