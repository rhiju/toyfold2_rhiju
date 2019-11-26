function [Emap,Jmap,theta] = scan_robot_arm_energy( sigma, theta0, dtheta )
%
% Make a heat map of 5-linker robot arm
%  with two 'driving' angles, and the other three
%  solved through kinematic closure
%
if ~exist( 'theta0','var') theta0 = 2*pi/5; end;
if ~exist( 'sigma','var') sigma = 1; end;
if ~exist( 'dtheta','var') dtheta = 0.1; end;
theta = [-pi:dtheta:pi];
Emap = nan * ones(length(theta),length(theta) );
Jmap = Emap;
tic
for i = 1:length( theta)
    for j = 1:length( theta)
        theta_driver = [theta(i), theta(j)];
        [val,theta_closed,J,  E_out, theta_out,J_out] = get_robot_arm_energy(theta_driver,theta0,sigma);
        % actually need to be careful to sum over *both* solutions
        Emap(i,j) = -log( sum( exp(-E_out)./J_out ) / sum(1./J_out) );
        Jmap(i,j) = 1./sum( 1./J_out );
    end
end
toc
% make plots
subplot(2,1,1);
imagesc( theta*180/pi, theta*180/pi, Emap )
axis image
colorbar( 'EastOutside' )
title( sprintf('Emap n_{link} = %d; theta=%8.3f; sigma=%8.3f',5,theta0,sigma) )
subplot(2,1,2);
imagesc( theta*180/pi, theta*180/pi, Jmap )
colorbar( 'EastOutside' )
axis image
set(gcf, 'PaperPositionMode','auto','color','white');
title( 'Jmap' );


