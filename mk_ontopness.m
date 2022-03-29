% see if there is ontopness in pattern

function [ontopness, count_ontop, ontop, below] = mk_ontopness(REQUIRED_FORM)

ontopness = false;
count_ontop = 0;

ontop = 0;
below = 0;

if ~isempty(REQUIRED_FORM)

    for i=1:(size(REQUIRED_FORM,1)-1)

        if any((REQUIRED_FORM(end-(i-1),:)-REQUIRED_FORM(end-i,:))~=0)

            index = find((REQUIRED_FORM(end-(i-1),:)-REQUIRED_FORM(end-i,:))~=0);

            if any((REQUIRED_FORM(end-(i-1),index).*REQUIRED_FORM(end-i,index))~=0)
                ontopness = true;
                count_ontop = count_ontop + 1;
                elements_form = unique(REQUIRED_FORM); elements_form(1) = [];
                [row1,~] = find(REQUIRED_FORM==elements_form(1));
                [row2,~] = find(REQUIRED_FORM==elements_form(2));
                if min(row1)<min(row2)
                    ontop  = elements_form(1);
                    below = elements_form(2);
                elseif min(row1)>min(row2)
                    ontop  = elements_form(2);
                    below = elements_form(1);
                end                
            end

        end
    end

end

end