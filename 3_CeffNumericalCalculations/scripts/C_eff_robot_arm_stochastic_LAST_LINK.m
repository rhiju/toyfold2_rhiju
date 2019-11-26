function C_eff = C_eff_robot_arm_stochastic_LAST_LINK( sigma, NITER, dr, dtheta, n_link, theta0 )
% Should be a little less choppy than naive RNAmake-style -- compute first
% three links of the robot arm, then compute probability density of
% last link

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;


x_all = zeros(2, NITER);
theta_all = zeros(1,NITER);

if ~exist( 'dr', 'var') dr = 0.01; end;
if ~exist( 'dtheta','var') dtheta = 0.01; end; % radians (cutoff for 'delta' function)
clf
tic
for i = 1:NITER
    theta = 0;
    x= [0,0];
    xpath = [0,0]';
    for n = 1:n_link-1
        theta = theta + theta0;
        theta = theta + sigma*randn(1); 
        x = x + [cos( theta) sin(theta) ];    
        xpath(:,n+1) = x;
    end
    x_all(:,i) = x;
    xpath(:,n_link+1)=[0,0];
    % unit vectors in coordinate frame at end of this third link
    vx = [cos(theta), sin(theta)];
    vy = [-sin(theta), cos(theta)];
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
    % how does origin coordinate frame map into this new frame?
    t_final = [0,0]-x; 
    r_final = norm(t_final);
    phi_final = atan2( (t_final*vy'), (t_final*vx') );
    theta_final = 2*pi - theta;
    C_eff_all(i) = 0;
    accepted(i) = 0;
    if abs(r_final-1)<(dr/2)
        d_angle = mod( (theta_final - phi_final) + pi, 2*pi ) - pi;
        if abs(d_angle) < dtheta/2   
            accepted(i) = 1;
            C_eff_all(i) = 2*pi *(1/dtheta)*(1/dr)*wrapped_gaussian(d_angle,sigma);
        end
    end    
end
C_eff = mean(C_eff_all);
toc

plot( x_all(1,:),x_all(2,:), 'b.' );
hold on
gp = find( accepted );
plot( x_all(1,gp),x_all(2,gp), 'r.' );
title( sprintf('NITER = %d; theta0 = %f; sigma = %f',NITER,theta0,sigma) )
make_polar_axes;