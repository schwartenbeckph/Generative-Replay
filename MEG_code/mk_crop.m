% crop non-zero parts out of columns/rows of all zeros
function [FORM_new, range_x, range_y] = mk_crop(FORM)

    for idx_form=1:size(FORM,3)

        min_x = min(find(sum(FORM(:,:,idx_form),1)~=0));
        max_x = max(find(sum(FORM(:,:,idx_form),1)~=0));

        min_y = min(find(sum(FORM(:,:,idx_form),2)~=0));
        max_y = max(find(sum(FORM(:,:,idx_form),2)~=0));
        
        if idx_form==1
            FORM_new = zeros(length(min_y:max_y),length(min_x:max_x),size(FORM,3));
        end

        FORM_new(:,:,idx_form) = FORM(min_y:max_y,min_x:max_x,idx_form);

        range_x(idx_form) = max_x - min_x + 1;
        range_y(idx_form) = max_y - min_y + 1;

    end

end