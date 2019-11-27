function [hist_forward,x_range,theta_range] = get_histogram(pts, L, dL, dtheta);
% [hist_forward,x_range,theta_range] = get_histogram(pts, L, dL, dtheta);

Nxy = 2*(L/dL)+1;
Ntheta = round((2*pi)/dtheta);
x_range = [0:(Nxy-1)]*dL-L;
theta_range = [0:(Ntheta-1)]*dtheta;
hist_forward = zeros( Nxy, Nxy, Ntheta );
for i = 1:size( pts, 1)
    hist_forward = update_histogram( pts(i,1:2), pts(i,3), hist_forward, L, dL, dtheta );
end
