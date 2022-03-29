% obtained matrices that reflect experienced transitions
% bricks_conn  Element connected to middle via besideness | middle Element | Element connected to middle via ontopness
% bricks_rel   Left element | ontop element | right element | below element
function [T_EE,T_ER,T_RE,T_EE_wrong,T_ER_wrong,T_RE_wrong,T_RR,T_RR_rev,T_RR_wrong] = mk_experienced_transitions(bricks_conn,bricks_rel,nstates,do_plot)

    if nargin<4
        do_plot = false;
    end
    
    %%%%% obtain experienced transitions between relations and elements %%%%%
    T_RE_temp  = zeros(nstates); % E1 | E2 | E3 | E4 | Left | Ontop | Right | Below
    T_ER_temp  = zeros(nstates);
    T_EE_temp  = zeros(nstates);
    
    T_ER = zeros(nstates);
    T_RE = zeros(nstates);
    
    T_RR       = zeros(nstates);
    T_RR_rev   = zeros(nstates);
    
    T_RE_wrong = zeros(nstates);
    T_ER_wrong = zeros(nstates);
    T_EE_wrong = zeros(nstates);
    T_RR_wrong = zeros(nstates);

    % make EE, ER and RE transitions:
    for idx_stim = 1:size(bricks_rel,1)

        if bricks_rel(idx_stim,1)~=0 && nstates>4
            T_RE_temp(5,bricks_rel(idx_stim,1)) = 1; % left to left brick     
            T_ER_temp(bricks_rel(idx_stim,1),7) = 1; % move right from left brick             
        end
        
        if bricks_rel(idx_stim,2)~=0 && nstates>4
            T_RE_temp(6,bricks_rel(idx_stim,2)) = 1; % ontop to ontop brick
            T_ER_temp(bricks_rel(idx_stim,2),8) = 1; % move below from ontop brick
        end
        
        if bricks_rel(idx_stim,3)~=0 && nstates>4
            T_RE_temp(7,bricks_rel(idx_stim,3)) = 1; % right to right brick
            T_ER_temp(bricks_rel(idx_stim,3),5) = 1; % move left from right brick            
        end
        
        if bricks_rel(idx_stim,4)~=0 && nstates>4
            T_RE_temp(8,bricks_rel(idx_stim,4)) = 1; % below to below brick
            T_ER_temp(bricks_rel(idx_stim,4),6) = 1; % move ontop from below brick
        end        
        
        if bricks_conn(idx_stim,1)~=0
            T_EE_temp(bricks_conn(idx_stim,2),bricks_conn(idx_stim,1)) = 1;
            T_EE_temp(bricks_conn(idx_stim,1),bricks_conn(idx_stim,2)) = 1;
        end
        
        if bricks_conn(idx_stim,3)~=0
            T_EE_temp(bricks_conn(idx_stim,2),bricks_conn(idx_stim,3)) = 1;
            T_EE_temp(bricks_conn(idx_stim,3),bricks_conn(idx_stim,2)) = 1;
        end        
        
    end
    
    if  nstates>4
        T_RE     = T_RE_temp;
        T_ER     = T_ER_temp;
    end
    T_EE     = T_EE_temp;    
    
    % make RR and RR back-and-forth (RR_rev) transitions:
    if size(bricks_rel,1)==1 && nstates>4
        
        if all(bricks_rel~=0)
            
            T_RR_rev(5,7) = 1; % left to right
            T_RR_rev(7,5) = 1; % right to left
            T_RR_rev(6,8) = 1; % ontop to below
            T_RR_rev(8,6) = 1; % below to ontop    

            if bricks_rel(1)==bricks_rel(2) % left element is ontop element => ontop to right and left to below
                T_RR(6,7) = 1; % ontop to right
                T_RR(5,8) = 1; % left to below
            elseif bricks_rel(2)==bricks_rel(3) % right element is ontop element => right to below or ontop to left
                T_RR(7,8) = 1; % right to below
                T_RR(6,5) = 1; % ontop to left
            elseif bricks_rel(3)==bricks_rel(4) % right element is below element => right to ontop or below to left
                T_RR(7,6) = 1; % right to ontop
                T_RR(8,5) = 1; % below to left
            elseif bricks_rel(1)==bricks_rel(4) % left element is below element => below to right or left to ontop
                T_RR(8,7) = 1; % below to right
                T_RR(5,6) = 1; % left to ontop
            end
            
        else
            
            if bricks_conn(1)==0 % ontopness relation probed
                T_RR_rev(6,8) = 1; % ontop to below
                T_RR_rev(8,6) = 1; % below to ontop                                
            elseif bricks_conn(3)==0 % besideness relation probed
                T_RR_rev(5,7) = 1; % left to right
                T_RR_rev(7,5) = 1; % right to left                
            end
        
        end
        
    end

    % make control sequences:
    if size(bricks_rel,1)==1
        
        missing_element = setdiff(1:4,bricks_rel);
        
        for idx_el=1:length(missing_element)
        
            if  nstates>4
                T_RE_wrong(5:end,missing_element(idx_el))           = 1;

                T_ER_wrong(missing_element(idx_el),5:end)           = 1;
            end

            T_EE_wrong(missing_element(idx_el),1:4)             = 1;
            T_EE_wrong(1:4,missing_element(idx_el))             = 1;
            T_EE_wrong(missing_element(idx_el),missing_element(idx_el)) = 0; 
        
        end
        
        if  nstates>4
            T_RR_wrong(5:end,5:end)                 = 1;
            T_RR_wrong(find(eye(size(T_RR_wrong)))) = 0; % don't include self-transitions, separate control
            T_RR_wrong                              = T_RR_wrong-(T_RR+T_RR_rev+T_RR');  
        end
        
        
    else    

        % make some other useful matrices
        if  nstates>4
            T_RE_wrong(5:end,1:4)                   = 1; 
            T_RE_wrong(find(eye(size(T_RE_wrong)))) = 0; % don't include self-transitions, separate control
            T_RE_wrong                              = T_RE_wrong-T_RE;

            T_ER_wrong(1:4,5:end)                   = 1;
            T_ER_wrong(find(eye(size(T_ER_wrong)))) = 0; % don't include self-transitions, separate control
            T_ER_wrong                              = T_ER_wrong-T_ER; 
        end

        T_EE_wrong(1:4,1:4)                     = 1;
        T_EE_wrong(find(eye(size(T_EE_wrong)))) = 0; % don't include self-transitions, separate control
        T_EE_wrong                              = T_EE_wrong-T_EE;            
    
    end
    
    if do_plot
        if  nstates>4
            figure,subplot(2,2,1),imagesc(T_RR),subplot(2,2,3),imagesc(T_RR_rev),subplot(2,2,4),imagesc(T_RR_wrong)            
            figure,subplot(2,2,1),imagesc(T_EE),subplot(2,2,3),imagesc(T_ER),subplot(2,2,4),imagesc(T_RE)
        end
        figure,subplot(2,2,1),imagesc(T_EE_wrong),subplot(2,2,3),imagesc(T_ER_wrong),subplot(2,2,4),imagesc(T_RE_wrong)

    end

end