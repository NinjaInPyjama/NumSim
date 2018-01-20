function v = onehot(ind, dim)
    % ONEHOT creates a (1 x dim) onehot vector with zeros everywhere except on index
    % ind   - index that is nonzero (is 1)
    % dim   - dimesion of the whole vertor
    
    if(size(ind,1) > size(ind,2))
        ind = ind';
    end
    v = zeros(dim,max(size(ind)));
    ind = ind + (0:size(ind,2)-1)*dim;
    v(ind) = 1;
    v= v';
end