function C_eff = C_eff_robot_arm_stochastic_GridOverlap( sigma, NITER, dL, dtheta, n_link, theta0, n_forward, MAKE_PLOT )
% C_eff = C_eff_robot_arm_stochastic_GridOverlap( sigma, NITER, dL, dtheta, n_link, theta0, n_forward, MAKE_PLOT )
% %
% Compute C_eff based on overlap of forward and reverse gridded histogram
%  for n-link polymer with unit length links
%
% INPUTS:
%  sigma  = width of wrapped Gaussian for angle (radians)
%  NITER  = number of independent samples
%  dL     = bin size in x,y  [Default = 0.1]
%  dtheta = bin size in theta (radians)
%  n_link = number of polymer links [Default 4]
%  theta0 = preferred angle at each link (radians) [Default = 2*pi/n_link]
%  n_forward = number of forward steps (n_reverse will be set to
%                   n_link - n_forward)  [Default = n_link/2]
%  MAKE_PLOT = make a plot of samples for visual feedback [Default = 0]
%
% OUTPUTS:
%  hist_reverse = 3D histogram in x, y, theta 
%  x_range      = bin_centers in (x, y) from -L to L with dL spacings
%  
% (C) R. Das, Stanford University 2019


if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;
if ~exist( 'n_forward','var') n_forward = round(n_link/2); end;
if ~exist( 'dL','var') dL = 0.1; end;
if ~exist( 'L','var') L = 4; end;
if ~exist( 'MAKE_PLOT','var') MAKE_PLOT = 0; end;

[hist_forward,x_range] = get_hist_forward(  NITER, sigma, theta0, dL, dtheta, L, n_forward );
n_reverse = n_link - n_forward;
[hist_reverse,x_range] = get_hist_reverse(  NITER, sigma, theta0, dL, dtheta, L, n_reverse );

C_eff = (2*pi) * sum(hist_forward(:)/NITER .* hist_reverse(:)/NITER/dtheta/dL^2);

if MAKE_PLOT
    subplot(2,1,1);
    imagesc( x_range,x_range,sum(hist_forward,3)' );
    colormap( 1 - copper(100) );
    make_polar_axes;
    
    subplot(2,1,2);
    imagesc( x_range,x_range,sum(hist_reverse,3)' );
    title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
    make_polar_axes;
end
