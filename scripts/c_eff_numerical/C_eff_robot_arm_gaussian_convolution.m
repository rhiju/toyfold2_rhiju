function C_eff_norm = C_eff_robot_arm_gaussian_convolution( n_link, theta0 )
% Essentially analytical -- propagate covariance matrix for each arm.
%  Makes assumption that sigma is small. 
% Note that output is C_eff* sigma^3.

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;

theta = 0;
x = [0,0];

% note funny convention -- origin is frame 1, first point is frame 2, etc. 
% used to match other computation scripts.
theta_path(1) = 0;
xpath = [0,0]';

for n = 1:n_link
    theta = theta + theta0;
    x = x + [cos( theta) sin(theta) ];
    theta_path(n+1) = theta;
    xpath(:,n+1) = x;
end
xpath

if theta0 == 2*pi/n_link; 
    assert( abs( principal_angle_radians(theta_path(end)) )<1e-5);
    assert( abs(xpath(1,end))<1e-5 );
    assert( abs(xpath(2,end))<1e-5 );
end;

cov_local = [0 0 0; 0 0 0; 0 0 1];
cov_total = zeros(3,3);
for n = 1:n_link
    % what is transform to get from local to global?
    [x,y,theta] = get_transform( xpath(:,end),theta_path(end),xpath(:,n),theta_path(n) );
    T = [cos(theta),sin(theta),-y*cos(theta)+x*sin(theta); -sin(theta), cos(theta), x*cos(theta)+y*sin(theta); 0,0,1 ];
    cov_global = T * cov_local * T';
    cov_total = cov_total + cov_global;
end
cov_total
C_eff_norm = (2*pi) * 1 / (2*pi)^(3/2) / sqrt(det( cov_total));

[x,y,theta] = get_transform( [0 0]', 0, xpath(:,end),theta_path(end));
C_eff_norm = C_eff_norm * exp( -0.5 * [x,y,theta] * inv(cov_total) * [x,y,theta]' );

function [x,y,theta] = get_transform( x_end,  theta_end, x_start, theta_start);
theta = principal_angle_radians( theta_end - theta_start );
v_x = [ cos(theta),-sin(theta)];
v_y = [ sin(theta), cos(theta)];
d_global = x_end-x_start;
x = v_x * d_global;
y = v_y * d_global;
