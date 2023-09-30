% normalise data
function X = mk_normalise(X)

%     figure,imagesc(X),colorbar

%     X = X./max(abs(X(:)));
    X = X./prctile(abs(X(:)),95);  
    
%     figure,imagesc(X),colorbar

end