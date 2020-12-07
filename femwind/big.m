function v=big(A)
% v=big(A)
% max abs element of nD matrix
    v=full(max(abs(A(:))));
end