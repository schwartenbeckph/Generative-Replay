% Recentre object in nxn-grid

function STIM = mk_recentre(FORM,put_where)

if nargin==1
    put_where = 'bottom_left'; % place in left bottom corner
end


STIM = zeros(size(FORM,1),size(FORM,2),size(FORM,3));

for i=1:size(FORM,3)
    
    n_grid = size(FORM(:,:,i),1); % matrix has to be square!

    stim = FORM(:,:,i);
    
    stim = stim(sum(stim,2)~=0,sum(stim,1)~=0); 
    
    if strcmp(put_where,'bottom_left')
        stim = [stim zeros(size(stim,1),n_grid-size(stim,2))];     
        stim = [zeros(n_grid-size(stim,1),size(stim,2)); stim];
    elseif strcmp(put_where,'top_left')
        stim = [stim zeros(size(stim,1),n_grid-size(stim,2))];     
        stim = [stim; zeros(n_grid-size(stim,1),size(stim,2))];   
    elseif strcmp(put_where,'middle')
        stim = [zeros(size(stim,1),ceil((n_grid-size(stim,2))/2)) stim zeros(size(stim,1),floor((n_grid-size(stim,2))/2))];     
        stim = [zeros(ceil((n_grid-size(stim,1))/2),size(stim,2)); stim; zeros(floor((n_grid-size(stim,1))/2),size(stim,2))]; 
    end
    
    STIM(:,:,i) = stim;
end


end