function [hist_forward,pts, x_range,theta_range] = get_hist_forward(  NITER, sigma, theta0, dL, dtheta, L, n_forward);
% [hist_forward,x_range,theta_range] = get_hist_forward(  NITER, sigma, theta0, dL, dtheta, L, n_forward );
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
%  n_forward = number of polymer links
%
% OUTPUTS:
%  hist_forward = 3D histogram in x, y, theta
%  x_range      = bin centers in  x (or y) from -L to L with dL spacings
%  theta_range  = bin centers in theta from 0 to 2 pi.%
%
% (C) R. Das, Stanford University 2019

tic
pts = get_pts_forward( NITER, sigma, theta0, n_forward );
[hist_forward,x_range,theta_range] = get_histogram(pts, L, dL, dtheta);
toc

