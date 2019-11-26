function val = get_signed_angle(a,b,c)
%
%     a
%      \  <----angle
%       b - c
%
% angle should come back from -pi to pi
%

v1 = a-b;
v2 = b-c;

% define a 'coordinate system' for v1
v1x = v1/norm(v1);
v1y = [-v1x(2),v1x(1)]';
val = atan2( v2'*v1y, v2'*v1x );
val = principal_angle_radians(val);
