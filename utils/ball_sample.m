function [X] = ball_sample(N, d)
%dropped coordinate method
X_aug = sphere_sample(N, d+2);
X = X_aug(:, 1:d);

end

