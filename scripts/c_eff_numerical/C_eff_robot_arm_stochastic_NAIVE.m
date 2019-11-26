function [C_eff,x_all,theta_all] = C_eff_robot_arm_stochastic_NAIVE( sigma, NITER, dr, dtheta, n_link, theta0 )
% This is RNAmake-style -- compute entire four-linker arm
% then see if end is close to origin and at right angle, within
% a tolerance that defines the 'area'.

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;

x_all = zeros(2, NITER);
theta_all = zeros(1,NITER);

if ~exist( 'dr', 'var') dr = 0.01; end;
if ~exist( 'dtheta','var') dtheta = 0.01; end; % radians (cutoff for 'delta' function)
tic
clf;
for i = 1:NITER
    theta = 0;
    x= [0,0];
    xpath = [0,0]';
    for n = 1:n_link
        theta = theta + theta0;
        theta = theta + sigma*randn(1); 
        x = x + [cos( theta) sin(theta) ];     
        xpath(:,n+1) = x;
    end
    x_all(:,i) = x;
    theta_all(:,i) = theta;
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
end
theta_all = mod( theta_all+pi, 2*pi ) - pi; % shift angle to be between -pi to pi
gp = find( (sum(x_all.^2) < dr.^2) & (abs(theta_all) < dtheta ) );
[length(gp) NITER]
std( principal_angle_radians( theta_all  ) )
sum( abs( principal_angle_radians( theta_all ) ) < dtheta )
sum((sum(x_all.^2) < dr.^2))
plot( x_all(1,:),x_all(2,:), 'b.' );
hold on
plot( x_all(1,gp),x_all(2,gp), 'r.' );
title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
C_eff = (2*pi) * (length(gp)/NITER) / (2*dtheta*pi*dr^2 );
toc

make_polar_axes;