function Esol = get_robot_arm_energy_from_theta_closed( theta_closed, theta0, sigma);
Esol = 0;
for j = 1:length(theta_closed )
    Esol = Esol - log_wrapped_gaussian(theta_closed(j)-theta0,sigma);
    if ( abs(theta_closed(j)-pi)<1e-5 ) Esol = Inf; end;
    if ( abs(theta_closed(j)+pi)<1e-5 ) Esol = Inf; end;
    if ( abs(theta_closed(j))<1e-5 )    Esol = Inf; end;
end
