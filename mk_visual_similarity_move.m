% assess visual similarity of stimuli in this task
% It's important that form1 and form2 are placed on left bottom corner of grid!

function vis_similarity_move = mk_visual_similarity_move(form1,form2)

 n_grid = 6;
 moves  = [+6,-1]'; % moves of translations, one to right (+6) or one up (-1)

 form1 = form1(form1~=0);
 form2 = form2(form2~=0);
 
 FORM1 = zeros(n_grid); FORM1(form1) = 1; FORM1 = mk_crop(FORM1);
 FORM2 = zeros(n_grid); FORM2(form2) = 1; FORM2 = mk_crop(FORM2);
 
 size1 = size(FORM1); size2 = size(FORM2);
 
 % find all translations of FORM1 in n_grid x n_grid Grid
 %  translations1 = n_grid-size1; % that's far too many. Rather, find translations where there's still overlap
 translations1     = size2-1; % find translations where there's still overlap
 max_translations1 = n_grid-size1;
 if translations1(1)>max_translations1(1)
     translations1(1)=max_translations1(1);
 end
  if translations1(2)>max_translations1(2)
     translations1(2)=max_translations1(2);
 end
 
 
%  [X1,Y1] = meshgrid(0:translations1(1),0:translations1(2));
 [X1,Y1] = meshgrid(0:translations1(2),0:translations1(1));
 TRANSLATIONS1 = cat(2,X1',Y1'); TRANSLATIONS1 = reshape(TRANSLATIONS1,[],2);
 
 TRANSLATIONS1 = TRANSLATIONS1*moves;
 
 all_forms1 = TRANSLATIONS1 + form1;
 
 % find all translations of FORM2 in n_grid x n_grid Grid
 %  translations2 = n_grid-size2; % that's far too many. Rather, find translations where there's still overlap
 translations2 = size1-1; % find translations where there's still overlap
 max_translations2 = n_grid-size2;
 if translations2(1)>max_translations2(1)
     translations2(1)=max_translations2(1);
 end
  if translations2(2)>max_translations2(2)
     translations2(2)=max_translations2(2);
 end
 
%  [X2,Y2] = meshgrid(0:translations2(1),0:translations2(2));
 [X2,Y2] = meshgrid(0:translations2(2),0:translations2(1));
 TRANSLATIONS2 = cat(2,X2',Y2'); TRANSLATIONS2 = reshape(TRANSLATIONS2,[],2);
 
 TRANSLATIONS2 = TRANSLATIONS2*moves;
 
 all_forms2 = TRANSLATIONS2 + form2;
 
 % now bring them all together, find overlap of all translated forms
 [a,b] = meshgrid(1:size(all_forms1,1),1:size(all_forms2,1));
 c = cat(2,a',b'); c = reshape(c,[],2);
 
%  test = [all_forms1(c(:,1),:) all_forms2(c(:,1),:)]
 
 vis_similarity_move = zeros(size(c,1),1);
 
 % a for loop here is not a very nice solution but works for now
 for i=1:size(c,1)
     % find proportion of overlap
     vis_similarity_move(i) = length(intersect(all_forms1(c(i,1),:),all_forms2(c(i,2),:)))/length(union(all_forms1(c(i,1),:),all_forms2(c(i,2),:)));
%      if vis_similarity_move(i)==1
%          error('stop')
%      end
 end
 
 vis_similarity_move = max(vis_similarity_move);

end