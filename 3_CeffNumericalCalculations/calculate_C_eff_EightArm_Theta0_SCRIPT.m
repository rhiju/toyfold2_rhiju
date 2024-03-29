n_link = 8;
theta0 = 0;
%% This is RNAmake-style.
NITER = 1000000;
sigma = [0.001 0.001 0.003 0.003 0.01 0.01 0.01 0.03 0.03 0.1 0.1 0.1 0.2 0.2 0.3 0.3 1];
C_eff_NAIVE = [];
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.01, 0.01, n_link, theta0 );
end
figure(2)
loglog( sigma, C_eff_NAIVE,'o');

% RNAmake style, but with tighter cutoffs.
C_eff_NAIVE2 = [];
NITER = 500000;
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE2(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.003, 0.003, n_link, theta0 );
end
figure(2)
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2],'o');

%% Mean C_eff at last link
C_eff_last_link = [];
NITER = 100000;
figure(1)
for i = 1:length(sigma)
    C_eff_last_link(i) = C_eff_robot_arm_stochastic_LAST_LINK( sigma(i), NITER, 0.01, 0.01, n_link, theta0);
end

%% Mean C_eff at last link
C_eff_last_link2 = [];
NITER = 100000;
figure(1)
for i = 1:length(sigma)
    C_eff_last_link2(i) = C_eff_robot_arm_stochastic_LAST_LINK( sigma(i), NITER, 0.003, 0.003, n_link, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o');


%% Mean C_eff after closure at two last links
C_eff_two_last_links = [];
NITER = 5000;
figure(1)
sigma_two_last_links = 10.^[-3:0.25:2];
for i = 1:length(sigma_two_last_links)
    C_eff_two_last_links(i) = C_eff_robot_arm_stochastic_TWO_LAST_LINKS( sigma_two_last_links(i), NITER, n_link, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on
plot( sigma_two_last_links, C_eff_two_last_links );


%% Harmonic expansion
C_eff_harmonic = [];
figure(1)
sigma_harmonic = 10.^[-3:0.25:2];
for i = 1:length(sigma_harmonic);
    C_eff_harmonic( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i), 300,10,n_link, theta0 ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, C_eff_harmonic);

%% Harmonic expansion
C_eff_harmonic2 = [];
figure(1)
for i = 1:length(sigma_harmonic);
    C_eff_harmonic2( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i),300,20,n_link, theta0 ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, [C_eff_harmonic; C_eff_harmonic2]);


%% Mining Minima 2 
C_eff_MM2 = [];
figure(1)
sigma_MM2 = 10.^[-3:0.1:2];
for i = 1:length(sigma_MM2);
    C_eff_MM2( i ) = C_eff_compute_MM2( sigma_MM2(i), n_link, theta0 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_MM2, C_eff_MM2);

%% Grid overlap
C_eff_grid_overlap1 = [];
C_eff_grid_overlap2 = [];
C_eff_grid_overlap3 = [];
C_eff_grid_overlap4 = [];
figure(1)
sigma_grid_overlap = 10.^[-3:0.1:2];
dL = 0.01; dtheta = 2*pi/100; theta0 = 0;
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

%% KDE-based overlap
C_eff_KDE_overlap1 = [];
C_eff_KDE_overlap2 = [];
C_eff_KDE_overlap3 = [];
C_eff_KDE_overlap4 = [];
figure(1); clf
sigma_KDE_overlap = 10.^[-3:0.1:2];
n_link = 8; theta0 = 0;
NITER = 1000;
for i = 1:length(sigma_KDE_overlap);
    C_eff_KDE_overlap1( i ) = C_eff_robot_arm_stochastic_KDE_Overlap( sigma_KDE_overlap(i), NITER, n_link, theta0, 1 );
    C_eff_KDE_overlap2( i ) = C_eff_robot_arm_stochastic_KDE_Overlap( sigma_KDE_overlap(i), NITER, n_link, theta0, 2 );
    C_eff_KDE_overlap3( i ) = C_eff_robot_arm_stochastic_KDE_Overlap( sigma_KDE_overlap(i), NITER, n_link, theta0, 3 );
    C_eff_KDE_overlap4( i ) = C_eff_robot_arm_stochastic_KDE_Overlap( sigma_KDE_overlap(i), NITER, n_link, theta0, 4 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_KDE_overlap, [C_eff_KDE_overlap1; C_eff_KDE_overlap2; C_eff_KDE_overlap3] );

%% In limit of large sigma, I have another way to compute, in terms of normal 2D Fourier expansion and Bessel functions...
% oh wait is this even bounded? If I increase x upper bound from 100 to
% 100000, I get higher and higher values.
dx = 0.001; x = [0:dx:1000];
C_eff_analytical = sum( besselj(0,x).^n_link .* x * dx ) / (2*pi);

%% Gaussian propagation
sigma_fine = 10.^[-3:0.1:2];
C_eff_gaussian_convolution = C_eff_robot_arm_gaussian_convolution(n_link, theta0) ./sigma_fine.^3; 

%%
set(figure(2),'position',[883 295 652 660]); clf;
h = plot( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
set(h,{'markerfacecolor'},get(h,'color'));
plot( sigma_two_last_links, C_eff_two_last_links, 'o','color',[0.3 0.7 0.3],'markerfacecolor',[0.3 0.7 0.3] );
plot( sigma_MM2, C_eff_MM2,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic2,'linew',1.5 );
plot( sigma_fine, C_eff_gaussian_convolution )
plot( sigma_fine, C_eff_analytical + 0*sigma_fine, 'color',[0.5 0.5 0.5] );
plot( sigma_grid_overlap, C_eff_grid_overlap1, 'o','color',[0.7 0.3 1]);
plot( sigma_grid_overlap, C_eff_grid_overlap2, 's','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap3, 'v','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap4, '^','color',[0.7 0.3 1]);
plot( sigma_KDE_overlap, C_eff_KDE_overlap1, 'o','color',[1 0.6 1]);
plot( sigma_KDE_overlap, C_eff_KDE_overlap2, 's','color',[1 0.6 1],'markerfacecolor',[1 0.6 1] );
plot( sigma_KDE_overlap, C_eff_KDE_overlap3, 'v','color',[1 0.6 1],'markerfacecolor',[1 0.6 1] );
plot( sigma_KDE_overlap, C_eff_KDE_overlap4, '^','color',[1 0.6 1]);

set(gca,'xscale','log','yscale','log');
legend( 'RNAmake-NAIVE','RNAmake-NAIVEfine','LastLink','LastLink-fine','TwoLastLinks',...
    'MiningMinima2','Harmonic to order 10 (Wang-Chirikjian)','Harmonic to order 20 (Wang-Chirikjian)','Gaussian conv', ...
    'Analytical (infinite sigma)','Grid Overlap1','Grid Overlap2','Grid Overlap3','Grid Overlap4',...
    'KDE Overlap1','KDE Overlap2','KDE Overlap3','KDE Overlap4');

xlabel( '\sigma (std. dev. of angle between arms)')
ylabel( 'C_{eff}' );
title( sprintf('C_{eff} predictions for effector with %d rigid arms to return to origin',n_link) );

set(gcf, 'PaperPositionMode','auto','color','white');
set( gca,'fontsize',14)

set(gcf, 'PaperPositionMode','auto','color','white');
set( gca,'fontsize',14);
ylim = get(gca,'ylim');
if ylim(1)<1e-12; set(gca,'ylim',[1e-12 ylim(2)] ); 
    set(h,'Location','SouthWest')
end;

