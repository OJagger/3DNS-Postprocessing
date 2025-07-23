function result = rowwise_slice(A, inds)
    nans = isnan(inds);
    inds(nans) = 1;
    inds = reshape(inds,1,[]);
    rows = (1:size(A,1));
    result = A(sub2ind(size(A), rows, inds));
    result(nans) = NaN;
end