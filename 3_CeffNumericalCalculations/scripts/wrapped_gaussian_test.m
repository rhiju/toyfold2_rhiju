%%
% Sanity check on wrapped Gaussian
sigma = 3;
dtheta = 0.005;
theta = -pi:dtheta:pi;
y = [];
for i = 1:length(theta); 
    y(i) = wrapped_gaussian(theta(i),sigma);
end
clf
plot( theta, y );
sum(y)*dtheta % should equal 1
