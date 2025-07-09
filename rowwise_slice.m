function result = rowwise_slice(A, inds)
    rows = (1:size(A,1))';
    result = A(sub2ind(size(A), rows, inds));
end