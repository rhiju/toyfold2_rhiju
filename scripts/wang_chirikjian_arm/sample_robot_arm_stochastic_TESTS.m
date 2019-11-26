%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model from Wang & Chirikjian, 3-jointed arms, fixed length 1, angle of 30
% degrees, sigma_theta = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
NITER = 1000; n_link = 3; theta0 = 30 * pi/180;
sigma = sqrt(0.1);
[x_all, theta_all] = sample_robot_arm_stochastic( n_link, theta0, sigma, NITER );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
NITER = 1000000; n_link = 4; theta0 = 90 * pi/180; sigma = 0.001;
[x_all, theta_all] = sample_robot_arm_stochastic( n_link, theta0, sigma, NITER );
cov( [x_all; theta_all]' )/(sigma^2) 
det( cov( [x_all; theta_all]' )/(sigma^2) )

figure(3);
NITER = 1000000; n_link = 4; theta0 = 90 * pi/180; sigma = 0.1;
[x_all, theta_all] = sample_robot_arm_stochastic( n_link, theta0, sigma, NITER );
cov( [x_all; theta_all]' )/(sigma^2) 
det( cov( [x_all; theta_all]' )/(sigma^2) )

figure(4);
NITER = 1000000; n_link = 4; theta0 = 90 * pi/180; sigma = sqrt(0.1);
[x_all, theta_all] = sample_robot_arm_stochastic( n_link, theta0, sigma, NITER );
cov( [x_all; theta_all]' )/(sigma^2) 
det( cov( [x_all; theta_all]' )/(sigma^2) )
