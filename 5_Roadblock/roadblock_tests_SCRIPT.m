% parameters
NITER = 5000;
sigma = 0.5;
theta0 = 0; %2*pi/10; 
n_max= 50; % number of links

%roadblock = [];  % no roadblock
roadblock = [5,0,  5]; % circular roadblock
trans_start = 0; rot_start = pi/2;
trans_final = [10,0]; rot_final = -pi/2;

%%
% main histogram loop --> forward/reverse
%pts_forward = {}; pts_reverse = {};
for n_link = 21:n_max
    fprintf( 'Doing %d of %d...\n',n_link,n_max);
    pts_forward{n_link} = get_pts_forward( NITER, sigma, theta0, n_link,trans_start,rot_start,roadblock);
    pts_reverse{n_link} = get_pts_reverse( NITER, sigma, theta0, n_link,trans_final,rot_final,roadblock );
end

%%
% plot histograms (x,y), marginalize over theta.
figure(1); 
dtheta = 2*pi/10; dL = 0.5; L = 25.0;
[hist_forward, hist_reverse] = get_all_histograms_from_pts( pts_forward, pts_reverse, dL, dtheta, L );

%%
% C_eff checks -- histogram/grid-based.
n_link = 20;
C_eff = [];
for n = 1:(n_link-1)
    C_eff(n) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
    fprintf( 'C_eff forward %d reverse %d [hist] ==> %f\n', n, n_link-n, C_eff(n) );
end

%%
all_C_eff = NaN * ones(n_link-1,n_link);
for n_link = 1:n_max
    for n = 1:(n_link-1)
        all_C_eff(n,n_link) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
        fprintf( 'C_eff forward %d reverse %d [hist] ==> %f\n', n, n_link-n, all_C_eff(n,n_link) );
    end
end
clf; plot( [1:n_max], all_C_eff','ko' );

%%
% C_eff checks -- KDE
n_link = 20;
C_eff1 = []; C_eff2 = [];
for n = 1:(n_link-1)
    C_eff1(n) = get_C_eff_from_pts( pts_forward{n}, pts_reverse{n_link-n} );
    fprintf( 'C_eff forward %d reverse %d [KDE1] ==> %f\n', n, n_link-n, C_eff1(n) );
    C_eff2(n) = get_C_eff_from_pts( pts_reverse{n_link-n}, pts_forward{n} );
    fprintf( 'C_eff forward %d reverse %d [KDE2] ==> %f\n', n, n_link-n, C_eff2(n) );
end
%%
% plot grid-based vs. C_eff
clf
plot( [1:n_link-1], [C_eff; C_eff1; C_eff2],'linew',1.5 );
title( sprintf('n_{link}=%d C_{eff}',n_link) );
legend( 'grid-based','KDE forward samples','KDE reverse samples' );
set(gcf, 'PaperPositionMode','auto','color','white');
xlabel( 'n_{forward}' ); ylabel( 'C_{eff}' );

%%
% backtrack to get samples.
figure(2); 
NSAMPLE = 10;
n_link = 20;  %n_link = 10;
plot_steps = 0;
use_mvksdensity = 1;
pts_reverse_already_transformed = 1;
all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans_final, rot_final, use_mvksdensity, plot_steps, trans_start, rot_start, pts_reverse_already_transformed, roadblock);

%%
% not updated yet to use trans_start, rot_start, and pre-transofrmed
% pts_reverse...
NSAMPLE = 20; all_trans = 3; trans_vec = [1,0];rot = 0;all_n_link=[5,10,20,30,40,50];
all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );

