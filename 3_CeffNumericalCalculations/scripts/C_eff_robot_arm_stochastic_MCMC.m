function [C_eff,x_all,theta_all] = C_eff_robot_arm_stochastic_MCMC( sigma, NITER, dr, dtheta, n_link, theta0, bias_strength )
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
theta_path = zeros(1,n_link)+theta0;
theta_path_prev = theta_path;
theta_step = 0.1;

e_prev = Inf;
clf; make_polar_axes;
counts = 0;
if ~exist( 'bias_strength' ) bias_strength = 1; end;
total_bias_weights = 0.0;
accepts = 0;
for i = 1:NITER
    theta_path = theta_path_prev;
    idx = randi(n_link);
    theta_path( idx ) = theta_path( idx ) + theta_step * randn(1);
    
    [x,theta,xpath] = get_path( theta_path );

    e = -sum( log_wrapped_gaussian( theta_path-theta0, sigma ) );
    bias_energy = 0.5 * bias_strength * norm(x);
    e = e + bias_energy;
    total_bias_weights = total_bias_weights + exp(bias_energy);
    if ( e > e_prev )
        if rand(1) > exp(e_prev-e) 
            theta_all(i)       = theta_all(i-1);
            x_all(:,i)         = x_all(:,i-1);
            weights_all(i) = weights_all(i-1);
            continue;
        end
    end
    accepts = accepts+1;

    theta_all(i) = theta;    
    x_all(:,i) = x';
    weights_all(i) = exp( bias_energy );
    if mod(i,100)==0;  plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); end;
    theta_path_prev = theta_path;    
    e_prev = e;
end
% theta_all = mod( theta_all+pi, 2*pi ) - pi; % shift angle to be between -pi to pi
% gp = find( (sum(x_all.^2) < dr.^2) & (abs(theta_all) < dtheta ) );
% plot( x_all(1,:),x_all(2,:), 'b.' );
% hold on
% plot( x_all(1,gp),x_all(2,gp), 'r.' );

plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] );
title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
fprintf( 'Accept rate: %f\n', accepts/NITER)
gp = find( (sum(x_all.^2) < dr.^2) & (abs(principal_angle_radians(theta_all)) < dtheta ) );
C_eff = (2*pi) * ( sum( weights_all(gp)) / sum( weights_all ) ) / (2*dtheta*pi*dr^2 );
fprintf( 'Number of hits: %d of %d. Bias strength %f. C_eff = %f\n', length(gp), NITER,bias_strength,C_eff);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end,theta_end,xpath] = get_path( theta_path );
xpath = [0,0]';
theta = 0;
x= [0,0];
theta = cumsum( theta_path );
xpath = [ [0,0]', cumsum( [cos(theta); sin(theta)],2 ) ];
x_end = xpath(:,end);
theta_end = theta(end);
