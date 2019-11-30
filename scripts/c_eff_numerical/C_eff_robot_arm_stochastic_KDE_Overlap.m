function C_eff = C_eff_robot_arm_stochastic_KDE_Overlap( sigma, NITER, n_link, theta0, n_forward, MAKE_PLOT )
% C_eff = C_eff_robot_arm_stochastic_KDE_Overlap( sigma, NITER, n_link, theta0, n_forward, MAKE_PLOT )
% %
% Compute C_eff based on overlap of forward and reverse KDEded histogram
%  for n-link polymer with unit length links
%
% INPUTS:
%  sigma  = width of wrapped Gaussian for angle (radians)
%  NITER  = number of independent samples
%  n_link = number of polymer links [Default 4]
%  theta0 = preferred angle at each link (radians) [Default = 2*pi/n_link]
%  n_forward = number of forward steps (n_reverse will be set to
%                   n_link - n_forward)  [Default = n_link/2]
%  MAKE_PLOT = make a plot of samples for visual feedback [Default = 0]
%
% OUTPUTS:
%  C_eff = Effective molarity to return to origin 
%  
% (C) R. Das, Stanford University 2019


if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;
if ~exist( 'n_forward','var') n_forward = round(n_link/2); end;
if ~exist( 'L','var') L = 4; end;
if ~exist( 'MAKE_PLOT','var') MAKE_PLOT = 0; end;

tic
pts_forward = get_pts_forward( NITER, sigma, theta0, n_forward );
n_reverse = n_link - n_forward;
pts_reverse = get_pts_reverse( NITER, sigma, theta0, n_reverse );
C_eff = get_C_eff_from_pts( pts_forward, pts_reverse );
toc
% 
% if MAKE_PLOT
%     subplot(2,1,1);
%     imagesc( x_range,x_range,sum(hist_forward,3)' );
%     colormap( 1 - copper(100) );
%     make_polar_axes;
%     
%     subplot(2,1,2);
%     imagesc( x_range,x_range,sum(hist_reverse,3)' );
%     title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
%     make_polar_axes;
% end
