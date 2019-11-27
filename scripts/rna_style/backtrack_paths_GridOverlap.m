% parameters
NITER = 20000;
sigma = 0.5;
theta0 = 0; %2*pi/10; 
n_max= 10; % number of links

%%
% main histogram loop --> forward
pts_forward = {}; pts_reverse = {};
for n_link = 1:n_max
    pts_forward{n_link} = get_pts_forward( NITER, sigma, theta0, n_link );
    pts_reverse{n_link} = get_pts_reverse( NITER, sigma, theta0, n_link );
end

%%
% plot histograms (x,y), marginalize over theta.
figure(1); 
dtheta = 2*pi/10; dL = 0.2; L = 12.0;
[hist_forward, hist_reverse] = get_all_histograms_from_pts( pts_forward, pts_reverse, dL, dtheta, L );

%%
% C_eff checks -- histogram
n_link = 10;
for n = 1:(n_link-1)
    C_eff(n) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
    fprintf( 'C_eff forward %d reverse %d ==> %f\n', n, n_link-n, C_eff(n) );
end

%%
% C_eff checks -- KDE
n_link = 10;
for n = 1:(n_link-1)
    C_eff(n) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
    fprintf( 'C_eff forward %d reverse %d ==> %f\n', n, n_link-n, C_eff(n) );
end

%%
% backtrack to get samples.
figure(2); 
NSAMPLE = 10;
n_link = 10;
use_mvksdensity = 1;
plot_steps = 0;
all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, use_mvksdensity, plot_stuff);


