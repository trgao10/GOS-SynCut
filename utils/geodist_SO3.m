function [ theta, K] = geodist_SO3( R1, R2 )

if abs(det(R1) - 1) > 1e-4  || abs(det(R2) - 1) > 1e-4
    disp('Input not SO(3).')
    disp( det(R2))
    theta = Inf;
    K = [];
    return
end
K = logm(real(R1)*real(R2)');
theta = norm( K, 'fro' )/sqrt(2) *sign(K(1,2));
K = K/theta;
