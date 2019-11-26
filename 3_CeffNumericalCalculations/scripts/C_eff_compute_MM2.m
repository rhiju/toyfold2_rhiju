function C_eff_total = C_eff_compute_MM2( sigma, n_link, theta0 );
% sigma = std dev. of wrapped gaussian
% n_link = number of links in robot arm

if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;
C_eff_total = 0;

% enumerate through solutions
for q = 1:(n_link-1)
    C_eff(q) = 0;
    theta_driver =  q* ( 2*pi/n_link ) * ones( 1, n_link-3 );
    
    [val,theta_closed,J] = get_robot_arm_energy(theta_driver,theta0,sigma);
    
    % need second derivatives (Hessian) -- just do it numerically.
    dtheta = 1e-5;
    for i = 1:(n_link-3)
        theta_di = theta_driver;
        theta_di(i) = theta_di(i) + dtheta;
        val_di = get_robot_arm_energy(theta_di,theta0,sigma);
        for j = 1:(n_link-3)
            theta_dj = theta_driver;
            theta_dj(j) = theta_dj(j) + dtheta;
            val_dj = get_robot_arm_energy(theta_dj,theta0,sigma);
            
            theta_di_dj = theta_driver;
            theta_di_dj(i) = theta_di_dj(i) + dtheta;
            theta_di_dj(j) = theta_di_dj(j) + dtheta;
            val_di_dj = get_robot_arm_energy(theta_di_dj,theta0,sigma);
            
            H(i,j) = ( val_di_dj - val_di - val_dj + val )/(dtheta*dtheta);
        end
    end
    if any( isnan( H(:) ) ); continue; end;
    if any( isinf( H(:) ) ); continue; end;
    %C_eff = (2*pi)*(2*pi)^((n_link-3)/2) / sqrt(det(H))  * exp(-val);
    [V,D] = eig(H);
    if any(D(:)<0); continue; end; % not a local minimum!
    C_eff(q) = (2*pi)*exp(-val)/J;
    for i = 1:(n_link-3)
        sigma_eff = 2/sqrt(D(i,i));
        C_eff(q) = C_eff(q) * sqrt(2*pi) * (sigma_eff/2) * erf( sqrt(2) * pi/sigma_eff);
    end
end

C_eff_total = sum( C_eff );
