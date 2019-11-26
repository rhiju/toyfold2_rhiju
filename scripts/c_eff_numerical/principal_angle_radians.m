function theta = principal_angle_radians( theta )
theta = mod( theta+pi, 2*pi ) - pi; % shift angle to be between -pi to pi
