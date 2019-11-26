function [val,theta_closed,J,E_out,theta_out,J_out] = get_robot_arm_energy(theta_driver,theta0,sigma)
% [val,theta_closed,J,E_out,theta_out,J_out] = get_robot_arm_energy(theta_driver,theta0,sigma)
%
%
E_out = [];
[theta_out,xpath_all,J_out] = find_closure_solutions( theta_driver );

if length( J_out ) == 0;
    theta_closed = [];
    J = NaN;
    val = Inf;
    return;
end

for q = 1:length( J_out )
    theta_closed = theta_out{q};
    E_out(q) = get_robot_arm_energy_from_theta_closed( theta_closed, theta0, sigma );
    [E_out(q), theta_closed];
end

[val,idx] = min( E_out );
theta_closed = theta_out{idx};
J = J_out(idx);
