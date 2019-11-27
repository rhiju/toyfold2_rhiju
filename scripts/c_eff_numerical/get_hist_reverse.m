function [hist_reverse,pts,x_range,theta_range] = get_hist_reverse(  NITER, sigma, theta0, dL, dtheta, L, n_reverse );
% [hist_reverse,x_range] = get_hist_reverse(  NITER, sigma, theta0, dL, dtheta, L, n_reverse );
%
% Compute gridded histogram
%  for n-link polymer
%
% INPUTS:
%  NITER  = number of independent samples
%  sigma  = width of wrapped Gaussian for angle (radians)
%  theta0 = preferred angle at each link (radians)
%  dL     = bin size in x,y
%  dtheta = bin size in theta (radians)
%  L      = max x/y
%  n_reverse = number of polymer links
%
% OUTPUTS:
%  hist_reverse = 3D histogram in x, y, theta
%  x_range      = bin centers in  x (or y) from -L to L with dL spacings
%  theta_range  = bin centers in theta from 0 to 2 pi.
%
% (C) R. Das, Stanford University 2019

tic
pts = get_pts_reverse( NITER, sigma, theta0, n_reverse );
[hist_reverse,x_range,theta_range] = get_histogram(pts, L, dL, dtheta);
toc

