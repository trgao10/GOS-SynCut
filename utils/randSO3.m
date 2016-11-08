function r = randSO3(sigma)

% sample random SO3 matrix with normal distribution on rotation angle theta

% sample angle
theta = rand*sigma;
while abs(theta) > pi
    theta = rand;
end

% sample rotation axis
k = rand(3,1);
k = k/norm(k);
% generate anti-symmetric matrix K
K = zeros(3);
K(1,2) = k(1);
K(1,3) = k(2);
K(2,3) = k(3);
K = K - K';

% take matrix exponential
r = expm(theta * K);
