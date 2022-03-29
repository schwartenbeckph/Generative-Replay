% see if there is besideness in pattern

function [besideness , count_beside, left, right] = mk_besideness(REQUIRED_FORM)

besideness = false;
count_beside = 0;

left  = 0;
right = 0;

if ~isempty(REQUIRED_FORM)
    
    REQUIRED_FORM = REQUIRED_FORM';

    for i=1:(size(REQUIRED_FORM,1)-1)

        if any((REQUIRED_FORM(end-(i-1),:)-REQUIRED_FORM(end-i,:))~=0)

            index = find((REQUIRED_FORM(end-(i-1),:)-REQUIRED_FORM(end-i,:))~=0);

            if any((REQUIRED_FORM(end-(i-1),index).*REQUIRED_FORM(end-i,index))~=0)
                besideness = true;
                count_beside = count_beside + 1;
                elements_form = unique(REQUIRED_FORM); elements_form(1) = [];
                [row1,~] = find(REQUIRED_FORM==elements_form(1));
                [row2,~] = find(REQUIRED_FORM==elements_form(2));
                if min(row1)<min(row2)
                    left  = elements_form(1);
                    right = elements_form(2);
                elseif min(row1)>min(row2)
                    left  = elements_form(2);
                    right = elements_form(1);
                end
            end

        end
    end

end

end