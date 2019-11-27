function [hist_forward,pts, x_range,theta_range] = get_hist_forward(  NITER, sigma, theta0, dL, dtheta, L, n_forward, use_mvksdensity );
% [hist_forward,x_range,theta_range] = get_hist_forward(  NITER, sigma, theta0, dL, dtheta, L, n_forward, use_mvksdensity );
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
if ~exist( 'use_mvksdensity', 'var' ); use_mvksdensity = 0; end;

tic
Nxy = 2*(L/dL)+1;
Ntheta = round((2*pi)/dtheta);
hist_forward = zeros( Nxy, Nxy, Ntheta );
x_range = [0:(Nxy-1)]*dL-L;
theta_range = [0:(Ntheta-1)]*dtheta;
pts = [];

for i = 1:NITER
    theta = 0;
    x = [0,0];
    xpath = [0,0]';
    for n = 1:n_forward
        theta = theta + theta0;
        theta = theta + sigma*randn(1);
        x = x + [cos( theta) sin(theta) ];
        xpath(:,n+1) = x;
    end
    pts = [pts; x,theta];
    hist_forward = update_histogram( x, theta, hist_forward, L, dL, dtheta );
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
end
toc
