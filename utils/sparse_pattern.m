function [Upsilon, Qzero] = sparse_pattern(S)
%SPARSE_PATTERN create the pattern for a matrix Q such that inv(Q) S has
%the same sparsity pattern as S
[n, m] = size(S');

[i, j] = find(S); % row index and column index of nonzero elements in a column vector
nb = length(i);
js = j + n*(0:nb-1)';
Sm = sparse(i, js, ones(size(js)), m, n*nb);


Q = reshape(1:(n^2), n, n);
Lam = ones(nb);

BQ = Sm*(kron(eye(nb), Q));
BLam = Sm*(kron(Lam, eye(n)));


%now find zero values
% [mi, mj, mv] = find((BLam==0) & (BQ~=0));
zero_ind = unique(BQ((BLam==0) & (BQ~=0)));

Qzero = Q;
Qzero(zero_ind) = 0;

Upsilon = ones(n);
Upsilon(zero_ind) = 0;

end

