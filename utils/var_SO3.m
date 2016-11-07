function v = var_SO3(C1, C2)

% input: C1, C2 are cell arrays containing SO3 matrices
% output: variation of the differences b/t C1 and C2 in euclidean space R^9

R = cellfun(@(x,y)real(x)*real(y)', C1, C2, 'UniformOutput', 0);
R = cat(4,R{:});
Rmean = mean(R, 4);
v = sum(reshape(R - repmat(Rmean, [1,1,1,size(R,4)]),[],1).^2);