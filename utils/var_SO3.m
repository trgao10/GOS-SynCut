function [v, Rmean] = var_SO3(C1, C2)

% input: C1, C2 are cell arrays containing SO3 matrices
% output: variation of the differences b/t C1 and C2 in euclidean space R^9

R = cellfun(@(x,y)real(x)*real(y)', C1, C2, 'UniformOutput', 0);
R = cat(3,R{:});
Rmean = mean(R, 3);
[U, ~, V] = svd(Rmean);
Rmean = U*V'; % this the minimizer of sum_i || r - r_i||_frob^2
v = sum(reshape(R - repmat(Rmean, [1,1,size(R,3)]),[],size(R,3)).^2,1);
v = sum(sqrt(v));