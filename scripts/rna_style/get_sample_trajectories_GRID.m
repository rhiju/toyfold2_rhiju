function all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
% all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );
% 
% Create grid of sample trajectories, moving target translation along
%  specified axes and ending at specified rotation -- can search
%  for pretty classes of trajectories (spirals, hearts, circles, etc.)
%
% INPUTS
%  NSAMPLE  = number of desired trajectories
%  all_n_link  = set of number of links in chain, sampled from top to
%                  bottom in grid display
%  all_trans   = set of translation magnitudes, sampled from left to right
%                  in grid display
%  trans_vec   = translation displacement direction (x,y)
%  rot      = target rotation (radians) [default 0]
%  theta0      = preferred angle at each link (radians)
%  sigma       = width of wrapped Gaussian for angle (radians)
%  pts_reverse = sampled (x,y,theta) for links built from origin in
%                    'reverse' direction
%  use_mvksdensity = use kernel density estimation [default: 1]
%  plot_steps = show intermediate steps & reverse histogram, pausing
%                       between each one [default 0]
%
% OUTPUT
%  all_xpaths = cell of [2 x n_link+1 x Nsample] generated paths shown as 
%                 grids. NaN marks paths that did not complete
%
% (C) R. Das, Stanford University, 2019

if ~exist( 'use_mvksdensity', 'var' ) use_mvksdensity = 1; end;
if ~exist( 'plot_steps', 'var') plot_steps = 0; end;

clf;
count = 0;
for i = 1:length(all_n_link)
    n_link = all_n_link(i);
    for j = 1:length(all_trans);
        count = count+1;
        trans = trans_vec*all_trans(j);
        subplot( length(all_n_link), length(all_trans), count );
        all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans,rot, use_mvksdensity, plot_steps);
        all_xpaths{i,j} = all_xpath;
        drawnow;
        set(gca,'ylim',[-9 9]);
        if ( j == 1 ); set(title( ['n_{link} = ',num2str(n_link)] ),'interpreter','tex'); end;
    end
end




