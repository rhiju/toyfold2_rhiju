function [theta_out,xpath_all,J_out] = find_closure_solutions( theta_driver )
% 
% theta_driver = 1, 2, ... N-3 input angles (in radians)
%
% Output
%  theta_out   = N-2, N-1, N angles needed to close N-link chain, in
%                     radians.
%                Note that multiple sets may be output (0 or 2)!
%  J_out       = det Jacobians associated with those closures |dq/dtheta_out|
%                     where q         = x_FINAL, y_FINAL, theta_FINAL
%                           theta_out = theta for N-1, N-1, N 

theta_out = {};
xpath_all = {};
J_out = [];

n_link = length( theta_driver )+3;
theta = 0;
x= [0,0];
xpath = [0,0]'; % start
L = 1;
for n = 1:(n_link-3)
    theta = theta + theta_driver(n);
    x = x + L * [cos( theta) sin(theta) ];
    xpath(:,n+1) = x;
end
% missing position of n_link-1. the rest are known, however.
xpath(:,n_link)  =[-L,0]; % will check this in a second.
xpath(:,n_link+1)=[0,0]; % end must return to start

d = xpath(:,n_link-2) - xpath(:,n_link);
D = norm(d);
% too long a distance to close chain?
if norm( D ) > (2*L); 
    return; % no solutions
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Need to define n-1, knowing positions
%       of n, n -2
%
%     n-1
%      [x]
%   L  /|\ L
%     / | \
%    x-----x <---  angle of triangle alpha
%    n  D  n-2
%

alpha = acos( 2/D );
% actually need to compute solution for alpha and -alpha.
bisect = (xpath(:,n_link-2) + xpath(:,n_link))/2;
d_vec = d/D; % direction vector
d_vec_perp = [d_vec(2), -d_vec(1)]; % direction perpendicular to d
d_perp = sqrt( L^2 - (D/2)^2);

% This new point defines three new bond angles  (n-3,n-2,n-1),
%   (n-2,n-1,n), and  (n-1, n, n_1).
% Compute jacobian involved in transforming these three angle
%  coordinates into x,y,theta coordinates of translation and rotation.
% This turns out to be 2 * area of the triangle defined by the new point
%  and its neighbors.
J =  d_perp * D;

count = 0;
for sgn = [-1 1] % flip triangle
    xpath(:,n_link-1) = bisect + sgn * d_perp * d_vec_perp';
    
    thetaA = get_signed_angle(xpath(:,n_link-3),xpath(:,n_link-2),xpath(:,n_link-1));
    thetaB = get_signed_angle(xpath(:,n_link-2),xpath(:,n_link-1),xpath(:,n_link  ));
    thetaC = get_signed_angle(xpath(:,n_link-1),xpath(:,n_link  ),xpath(:,n_link+1));
    
    count = count+1;
    theta_out{count} = [theta_driver, thetaA, thetaB, thetaC];
    xpath_all{count} = xpath;
    J_out(count) = J;
end



