function C_eff = C_eff_robot_arm_stochastic_TWO_LAST_LINKS( sigma, NITER, n_link, theta0 )
% Should be a even less choppy than naive RNAmake-style -- compute first
% few links of the robot arm, then compute geometry of last two links
% needed to complete chain; tally up C_eff, taking careful account of
% Jacobian.

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;

x_all = zeros(2, NITER);
theta_all = zeros(1,NITER);

clf
tic
C_eff_all = zeros(1,NITER);
L = 1;
theta_min = 1e-3;
for i = 1:NITER
    for n = 1:(n_link-3);        theta_driver(n) = theta0 + sigma*randn(1);     end
    
    [theta_out,xpath_all,J_out] = find_closure_solutions( theta_driver );

    for q = 1:length( J_out ) 
        theta_closed = theta_out{q}(end-2:end);
        if ( abs(theta_closed(1)) < theta_min ) continue; end;
        if ( abs(theta_closed(2)) < theta_min ) continue; end;
        if ( abs(theta_closed(3)) < theta_min ) continue; end;
        xpath = xpath_all{q};
        J = J_out(q);
        if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',1,'color',[0.5,0.5,0.5] ); hold on; end;
        C_eff_all(i) = C_eff_all(i) +  2*pi  ...
            * wrapped_gaussian(theta_closed(1)-theta0,sigma) ...
            * wrapped_gaussian(theta_closed(2)-theta0,sigma) ...
            * wrapped_gaussian(theta_closed(3)-theta0,sigma) ...
            / J;
    end
end
C_eff = mean(C_eff_all);
toc

title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
make_polar_axes;


