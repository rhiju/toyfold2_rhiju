% parameters
NITER = 5000;
sigma = 0.5;
theta0 = 0; %2*pi/10; 
n_max= 50; % number of links

%%
% main histogram loop --> forward
pts_forward = {}; pts_reverse = {};
for n_link = 1:n_max
    fprintf( 'Doing %d of %d...\n',n_link,n_max);
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
n_link = 50;
for n = 1:(n_link-1)
    C_eff(n) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
    fprintf( 'C_eff forward %d reverse %d [hist] ==> %f\n', n, n_link-n, C_eff(n) );
end

%%
% C_eff checks -- KDE
n_link = 50;
for n = 1:(n_link-1)
    C_eff1(n) = get_C_eff_from_pts( pts_forward{n}, pts_reverse{n_link-n} );
    fprintf( 'C_eff forward %d reverse %d [KDE1] ==> %f\n', n, n_link-n, C_eff1(n) );
    C_eff2(n) = get_C_eff_from_pts( pts_reverse{n_link-n}, pts_forward{n} );
    fprintf( 'C_eff forward %d reverse %d [KDE2] ==> %f\n', n, n_link-n, C_eff2(n) );
end

%%
% backtrack to get samples.
figure(2); 
NSAMPLE = 2;
n_link = 50;  %n_link = 10;
trans = [0,0];
rot = 0;
plot_steps = 0;
use_mvksdensity = 1;
all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans, rot, use_mvksdensity, plot_steps);

%%
% for fun.
clf;
n_link = 8;
NSAMPLE = 100;
trans = [3,0];
rot = 0;
plot_steps = 0;use_mvksdensity = 1;
all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans,rot, use_mvksdensity, plot_steps);

%%
NSAMPLE = 20; all_trans = 3; trans_vec = [1,0];rot = 0;all_n_link=[5,10,20,30,40,50];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );

%%
% do a big plot
% if ~exist( 'Figures','dir'); mkdir( 'Figures' ); end;
NSAMPLE = 10;
all_n_link = [4:2:10];
all_trans  = [-8 -4 -2 0 2 4 8]; 

rot = 0;
trans_vec  = [0 1];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
export_fig( 'Figures/lateral_GRID.pdf' );


rot = 0;
trans_vec = [1 0];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
export_fig( 'Figures/axis_GRID.pdf' );


rot = pi;
trans_vec  = [0 1];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
export_fig( 'Figures/lateral_reversedir_GRID.pdf' );


rot = pi;
trans_vec  = [1 0];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
export_fig( 'Figures/axis_reversedir_GRID.pdf' );


