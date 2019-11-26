function C_eff = C_eff_compute_semianalytic( sigma, theta0, dtheta);
if ~exist( 'dtheta','var') dtheta = 0.1; end;
[Emap,Jmap,theta] = scan_robot_arm_energy( sigma, theta0, dtheta );

C = exp(-Emap)./Jmap;
C( find( isnan(C(:)) ) ) = 0;
C( find( isinf(C(:)) ) ) = 0;
C_eff = 2*pi*sum(sum(C)) * dtheta.^2;



