function theta_closed = find_minima( sigma, NITER, n_link, theta0 )
%  theta_closed = find_minima( sigma, NITER, n_link )
%
% Should be a even less choppy than naive RNAmake-style -- compute first
% few links of the robot arm, then compute geometry of last two links
% needed to complete chain; tally up C_eff, taking careful account of
% Jacobian.

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;

for n = 1:NITER
    p = randn(1, n_link-3 );
    theta_driver = fminsearch( @get_robot_arm_energy, p, [], theta0,sigma );
    [val,theta_closed] = get_robot_arm_energy(theta_driver,theta0,sigma);
    %fprintf( 'theta_closed %f\n', theta_closed )
    %val   
    principal_angle_radians(theta_closed)*180/pi
end


%title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
%make_polar_axes;


