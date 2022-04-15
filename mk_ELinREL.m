
function REL_RSM = mk_ELinREL(stim_Indices,ElR_stim)

REL_RSM = zeros(length(stim_Indices),length(stim_Indices));
    
for idx_col=1:length(stim_Indices)
   for idx_row=1:length(stim_Indices)
       REL_RSM(idx_row,idx_col) = mean(ElR_stim(stim_Indices(idx_row),:) == ElR_stim(stim_Indices(idx_col),:));
   end
end
    
end