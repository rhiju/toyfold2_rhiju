%% Force preferred angle to be theta0 = 0 -- so default geometry is a straight rod!
theta0 = 0;

%% This is RNAmake-style.
NITER = 1000000;
sigma = [0.001 0.001 0.003 0.003 0.01 0.01 0.01 0.03 0.03 0.1 0.1 0.1 0.2 0.2 0.3 0.3 1];
C_eff_NAIVE = [];
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.01, 0.01, 5, theta0 );
end
figure(2)
loglog( sigma, C_eff_NAIVE,'o');

% RNAmake style, but with tighter cutoffs.
C_eff_NAIVE2 = [];
NITER = 500000;
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE2(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.003, 0.003, 5, theta0 );
end
figure(2)
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2],'o');

% Mean C_eff at last link
C_eff_last_link = [];
NITER = 100000;
figure(1)
for i = 1:length(sigma)
    C_eff_last_link(i) = C_eff_robot_arm_stochastic_LAST_LINK( sigma(i), NITER, 0.01, 0.01, 5, theta0);
end

% Mean C_eff at last link
C_eff_last_link2 = [];
NITER = 100000;
figure(1)
for i = 1:length(sigma)
    C_eff_last_link2(i) = C_eff_robot_arm_stochastic_LAST_LINK( sigma(i), NITER, 0.003, 0.003, 5, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o');


% Mean C_eff after closure at two last links
C_eff_two_last_links = [];
NITER = 50000;
figure(1)
sigma_two_last_links = 10.^[-3:0.1:2];
for i = 1:length(sigma_two_last_links)
    C_eff_two_last_links(i) = C_eff_robot_arm_stochastic_TWO_LAST_LINKS( sigma_two_last_links(i), NITER, 5, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on
plot( sigma_two_last_links, C_eff_two_last_links );


% Harmonic expansion
C_eff_harmonic = [];
figure(1)
sigma_harmonic = 10.^[-3:0.25:2];
for i = 1:length(sigma_harmonic);
    C_eff_harmonic( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i), 300,10,5,theta0 ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, C_eff_harmonic);

% Harmonic expansion
C_eff_harmonic2 = C_eff_harmonic;
figure(1)
for i = 1:length(sigma_harmonic);
    C_eff_harmonic2( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i),300,20,5, theta0 ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, [C_eff_harmonic; C_eff_harmonic2]);

% Mining Minima 2 
C_eff_MM2 = [];
figure(1)
sigma_MM2 = 10.^[-3:0.05:2];
for i = 1:length(sigma_MM2);
    C_eff_MM2( i ) = C_eff_compute_MM2( sigma_MM2(i), 5, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_MM2, C_eff_MM2);

% Semianalytic (2D integral over driver angles)
C_eff_semianalytic = [];
figure(1)
sigma_semianalytic = 10.^[-1.5:0.1:2];
dtheta=0.05;
for i = 1:length(sigma_semianalytic);
    C_eff_semianalytic( i ) = C_eff_compute_semianalytic( sigma_semianalytic(i), theta0, dtheta );
end
figure(2); clf
loglog( sigma_two_last_links, C_eff_two_last_links, 'o'); hold on;
plot( sigma_semianalytic, C_eff_semianalytic);

%% Grid overlap
C_eff_grid_overlap1 = [];
C_eff_grid_overlap2 = [];
C_eff_grid_overlap3 = [];
C_eff_grid_overlap4 = [];
figure(1)
sigma_grid_overlap = 10.^[-3:0.1:2];
dL = 0.01; dtheta = 2*pi/100;n_link = 5; theta0 = 0;
NITER = 10000;
for i = 1:length(sigma_grid_overlap);
    C_eff_grid_overlap1( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 1 );
    C_eff_grid_overlap2( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 2 );
    C_eff_grid_overlap3( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 3 );
    C_eff_grid_overlap4( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 4 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_grid_overlap, [C_eff_grid_overlap1; C_eff_grid_overlap2; C_eff_grid_overlap3] );

%% In limit of large sigma, I have another way to compute, in terms of normal 2D Fourier expansion and Bessel functions...
dx = 0.001; x = [0:dx:1000];
C_eff_analytical = sum( besselj(0,x).^5 .* x * dx ) / (2*pi);

%% Gaussian propagation
C_eff_gaussian_convolution = C_eff_robot_arm_gaussian_convolution(5,theta0) ./sigma_fine.^3; 

%% make final plot
set(figure(2),'position',[883 295 652 660]); clf;
h = plot( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
set(h,{'markerfacecolor'},get(h,'color'));
plot( sigma_two_last_links, C_eff_two_last_links, 'o','color',[0.3 0.7 0.3],'markerfacecolor',[0.3 0.7 0.3] );
plot( sigma_MM2, C_eff_MM2,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic2,'linew',1.5 );
plot( sigma_fine, C_eff_gaussian_convolution )
plot( sigma_fine, C_eff_analytical + 0*sigma_fine, 'color',[0.5 0.5 0.5] );
plot( sigma_semianalytic, C_eff_semianalytic, 'color','k' );
plot( sigma_grid_overlap, C_eff_grid_overlap1, 'o','color',[0.7 0.3 1]);
plot( sigma_grid_overlap, C_eff_grid_overlap2, 's','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap3, 'v','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap4, '^','color',[0.7 0.3 1]);


set(gca,'xscale','log','yscale','log');
h = legend( 'RNAmake-NAIVE','RNAmake-NAIVEfine','LastLink','LastLink-fine','TwoLastLinks',...
    'MiningMinima2','Harmonic to order 10 (Wang-Chirikjian)','Harmonic to order 20 (Wang-Chirikjian)','Gaussian conv', ...
    'Analytical (infinite sigma)','Semianalytical (2D integral)','Grid Overlap1','Grid Overlap2','Grid Overlap3','Grid Overlap4');
xlabel( '\sigma (std. dev. of angle between arms)')
ylabel( 'C_{eff}' );
title( 'C_{eff} predictions for effector with 5 rigid arms to return to origin [theta_0 = 0]' );

set(gcf, 'PaperPositionMode','auto','color','white');
set( gca,'fontsize',14);
ylim = get(gca,'ylim');
if ylim(1)<1e-10; set(gca,'ylim',[1e-10 ylim(2)] ); 
    set(h,'Location','SouthWest')
end;
