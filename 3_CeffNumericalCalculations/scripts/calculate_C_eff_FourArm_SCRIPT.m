%% This is RNAmake-style.
NITER = 1000000;
sigma = [0.001 0.001 0.003 0.003 0.01 0.01 0.01 0.03 0.03 0.1 0.1 0.1 0.2 0.2 0.3 0.3 1];
C_eff_NAIVE = [];
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.01, 0.01 );
end
figure(2)
loglog( sigma, C_eff_NAIVE,'o');

%% RNAmake style, but with tighter cutoffs.
C_eff_NAIVE2 = [];
NITER = 5000000;
figure(1)
for i = 1:length(sigma)
    C_eff_NAIVE2(i) = C_eff_robot_arm_stochastic_NAIVE( sigma(i), NITER, 0.003, 0.003 );
end
figure(2)
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2],'o');

%% Mean C_eff at last link
C_eff_last_link = [];
NITER = 10000;
figure(1)
for i = 1:length(sigma)
    C_eff_last_link(i) = C_eff_robot_arm_stochastic_LAST_LINK( sigma(i), NITER, 0.01, 0.01 );
end

%% Mean C_eff at two last links
C_eff_two_last_links = [];
NITER = 5000;
figure(1)
sigma_two_last_links = 10.^[-3:0.25:2];
for i = 1:length(sigma_two_last_links)
    C_eff_two_last_links(i) = C_eff_robot_arm_stochastic_TWO_LAST_LINKS( sigma_two_last_links(i), NITER, 4 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on
plot( sigma_two_last_links, C_eff_two_last_links );


%% Mean C_eff at last link
C_eff_two_last_links = [];
NITER = 1000;
figure(1)
for i = 1:length(sigma)
    C_eff_two_last_links(i) = C_eff_robot_arm_stochastic_TWO_LAST_LINKS( sigma(i), NITER );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2; C_eff_two_last_links],'o');


%% Semianalytical integration
sigma_fine = 10.^[-3:0.1:2];
C_eff_semianalytical = [];
for i = 1:length(sigma_fine)
    C_eff_semianalytical(i) = C_eff_robot_arm_semianalytical( sigma_fine(i) );
end
C_eff_semianalytical2 = [];
for i = 1:length(sigma_fine)
    C_eff_semianalytical2(i) = C_eff_robot_arm_semianalytical( sigma_fine(i), 1e-5 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on
plot( sigma_fine, [C_eff_semianalytical; C_eff_semianalytical2]);


%% Harmonic expansion
C_eff_harmonic = [];
figure(1)
sigma_harmonic = 10.^[-3:0.25:2];
for i = 1:length(sigma_harmonic);
    C_eff_harmonic( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i) ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, C_eff_harmonic);


%% Harmonic expansion
C_eff_harmonic2 = [];
figure(1)
for i = 1:length(sigma_harmonic);
    C_eff_harmonic2( i ) = real( C_eff_SE2_harmonic_expansion( sigma_harmonic(i),300,20 ) );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_harmonic, C_eff_harmonic2);


%% Mining Minima 2
C_eff_MM2 = [];
figure(1)
sigma_MM2 = 10.^[-3:0.1:2];
for i = 1:length(sigma_MM2);
    C_eff_MM2( i ) = C_eff_compute_MM2( sigma_MM2(i), 4 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_MM2, C_eff_MM2);


%% Grid overlap
C_eff_grid_overlap1 = [];
C_eff_grid_overlap2 = [];
C_eff_grid_overlap3 = [];
figure(1)
sigma_grid_overlap = 10.^[-3:0.1:2];
dL = 0.01; dtheta = 2*pi/100;n_link = 4; theta0 = 2*pi/n_link;
NITER = 10000;
for i = 1:length(sigma_grid_overlap);
    C_eff_grid_overlap1( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 1 );
    C_eff_grid_overlap2( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 2 );
    C_eff_grid_overlap3( i ) = C_eff_robot_arm_stochastic_GridOverlap( sigma_grid_overlap(i), NITER, dL, dtheta, n_link, theta0, 3 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_grid_overlap, [C_eff_grid_overlap1; C_eff_grid_overlap2; C_eff_grid_overlap3] );


%% Markov Chain Monte Carlo (with 'L1' bias)
C_eff_MCMC1 = [];
C_eff_MCMC2 = [];
C_eff_MCMC3 = [];
figure(1)
dr = 0.01; dtheta = 0.01;n_link = 4; theta0 = 2*pi/n_link;
NITER = 10000;
sigma_MCMC =  [0.001 0.001 0.003 0.003 0.01 0.01 0.01 0.03 0.03 0.1 0.1 0.1 0.2 0.2 0.3 0.3 0.5 0.5 1 1 2 2 3 3 5 5];
for i = 1:length(sigma_MCMC);
    C_eff_MCMC1( i ) = C_eff_robot_arm_stochastic_MCMC( sigma_MCMC(i), NITER, dr, dtheta, n_link, theta0, 0 );
    C_eff_MCMC2( i ) = C_eff_robot_arm_stochastic_MCMC( sigma_MCMC(i), NITER, dr, dtheta, n_link, theta0, 10 );
    C_eff_MCMC3( i ) = C_eff_robot_arm_stochastic_MCMC( sigma_MCMC(i), NITER, dr, dtheta, n_link, theta0, 100 );
end
figure(2); clf
loglog( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
plot( sigma_MCMC, [C_eff_MCMC1; C_eff_MCMC2; C_eff_MCMC3] );

%% In limit of large sigma, I have another way to compute, in terms of normal 2D Fourier expansion and Bessel functions...
% oh wait is this even bounded? If I increase x upper bound from 100 to
% 100000, I get higher and higher values.
dx = 0.001; x = [0:dx:100];
C_eff_analytical = sum( besselj(0,x).^4 .* x * dx ) / (2*pi);
dx = 0.001; x = [0:dx:100000];
C_eff_analytical2 = sum( besselj(0,x).^4 .* x * dx ) / (2*pi);
clear x

%% Gaussian propagation
C_eff_gaussian_convolution = 1/sqrt(8*pi) * 1./sigma_fine.^3;
hold on

%% Make plots
set(figure(2),'position',[883 295 652 660]); clf;
loglog( sigma_fine, C_eff_semianalytical, 'k', 'linew',1 ); hold on
loglog( sigma_fine, C_eff_semianalytical2, 'k', 'linew',1.5 ); hold on
h = plot( sigma, [C_eff_NAIVE; C_eff_NAIVE2; C_eff_last_link; C_eff_last_link2],'o'); hold on;
set(h,{'markerfacecolor'},get(h,'color'));
plot( sigma_two_last_links, C_eff_two_last_links, 'o','color',[0.3 0.7 0.3],'markerfacecolor',[0.3 0.7 0.3] );
plot( sigma_MM2, C_eff_MM2,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic,'linew',1.5 );
plot( sigma_harmonic, C_eff_harmonic2,'linew',1.5 );
plot( sigma_fine, C_eff_gaussian_convolution )
plot( sigma_fine, C_eff_analytical + 0*sigma_fine, 'color',[0.5 0.5 0.5] );
plot( sigma_fine, C_eff_analytical2 + 0*sigma_fine, 'linew',1.5,'color',[0.5 0.5 0.5] );
plot( sigma_grid_overlap, C_eff_grid_overlap1, 'o','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap2, 's','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_grid_overlap, C_eff_grid_overlap3, 'v','color',[0.7 0.3 1],'markerfacecolor',[0.7 0.3 1] );
plot( sigma_MCMC, C_eff_MCMC1, 'o','color',[1.0 0.5 0.4]);
plot( sigma_MCMC, C_eff_MCMC2, 's','color',[1.0 0.5 0.4]);
plot( sigma_MCMC, C_eff_MCMC3, 'v','color',[1.0 0.5 0.4]);
set(gca,'xscale','log','yscale','log');
legend( 'Semianalytical','Semianalytical-fine','RNAmake-NAIVE','RNAmake-NAIVEfine','LastLink','LastLink-fine','TwoLastLinks',...
    'MiningMinima2','Harmonic to order 10 (Wang-Chirikjian)','Harmonic to order 20 (Wang-Chirikjian)','Gaussian conv', ...
    'Analytical (infinite sigma)','Analytical (infinite sigma) -fine','Grid Overlap1','Grid Overlap2','Grid Overlap3','MCMC nobias','MCMC bias 10','MCMC bias 100');
xlabel( '\sigma (std. dev. of angle between arms, around 90\circ)')
ylabel( 'C_{eff}' );
title( 'C_{eff} predictions for effector with 4 rigid arms to return to origin' );

set(gcf, 'PaperPositionMode','auto','color','white');
set( gca,'fontsize',14)


%% Output for comparison to other ToyFold2 enthusiasts
sigma_test = [1.000e-03 3.594e-03 1.292e-02 4.642e-02 1.668e-01 5.995e-01 2.154e+00 7.743e+00 2.783e+01 1.000e+02];
titles = {'Semianalytical','Semianalytical-fine','RNAmake-NAIVE','RNAmake-NAIVEfine','LastLink','LastLink-fine','TwoLastLinks',...
    'MiningMinima2','SE(2) Harmonic to order 10 (Wang-Chirikjian)','SE(2) Harmonic to order 20 (Wang-Chirikjian)','Gaussian conv', ...
    'Analytical (infinite sigma)','Analytical (infinite sigma) -fine','Grid Overlap1','Grid Overlap2','Grid Overlap3'};
sigma_val = {sigma_fine,sigma_fine,sigma,sigma,sigma,sigma,sigma_two_last_links,...
    sigma_MM2,sigma_harmonic,sigma_harmonic,sigma_fine... 
    sigma_fine,sigma_fine,sigma_grid_overlap,sigma_grid_overlap,sigma_grid_overlap};
Ceff_val = {C_eff_semianalytical,C_eff_semianalytical2,C_eff_NAIVE, C_eff_NAIVE2, C_eff_last_link, C_eff_last_link2,C_eff_two_last_links,...
    C_eff_MM2,C_eff_harmonic,C_eff_harmonic2,C_eff_gaussian_convolution,...
    C_eff_analytical + 0*sigma_fine,C_eff_analytical2 + 0*sigma_fine,C_eff_grid_overlap1,C_eff_grid_overlap2,C_eff_grid_overlap1};
output_comparison_table( sigma_test, titles, sigma_val, Ceff_val );
