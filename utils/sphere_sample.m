function [X] = sphere_sample(N, d)
%dropped coordinate method
X = randn(N, d);
normU = sqrt(sum(X.^2, 2));
X = X ./ normU;

end

