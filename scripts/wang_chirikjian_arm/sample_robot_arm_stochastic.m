function [x_all, theta_all] = sample_robot_arm_stochastic( n_link, theta0, sigma, NITER );
% model from Wang & Chirikjian, 3-jointed arms, fixed length 1, angle of 30
% degrees, sigma_theta = 0.1
x_all = zeros(2, NITER);
theta_all = zeros(1,NITER);
for i = 1:NITER
    theta = 0;
    x= [0,0];
    for n = 1:n_link
        theta = theta + theta0;
        theta = theta + sigma*randn(1); 
        x = x + [cos( theta) sin(theta) ];      
    end
    x_all(:,i) = x;
    theta_all(:,i) = theta;
end
clf;
plot( x_all(1,:),x_all(2,:), '.' );
hold on
title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )

make_polar_axes;