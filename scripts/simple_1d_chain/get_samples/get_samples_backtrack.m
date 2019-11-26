NITER = 100000;
n_link = 6;
x = [0, 0.1*randn( 1, n_link-1 )]; % deviation of each link position from its ideal value.
sep_ideal = [-1,1,-1,1,-1,1];
x_ideal = [0,1,0,1,0,1];
% metropolis chain monte carlo
x_all = zeros(NITER,n_link);
for i = 1:NITER    
    x(1) = 0;
    x(2) = sample_based_on_neighbors( x(1), 1, x(1), sqrt(5) );
    x(3) = sample_based_on_neighbors( x(2), 1, x(1), sqrt(4) );
    x(4) = sample_based_on_neighbors( x(3), 1, x(1), sqrt(3) );
    x(5) = sample_based_on_neighbors( x(4), 1, x(1), sqrt(2) );
    x(6) = sample_based_on_neighbors( x(5), 1, x(1), sqrt(1) );    
    x = x + x_ideal;
    x_all(i,:) = x;    
end

subplot(2,1,1);
plot(x_all);
xlabel( 'Sample' ); ylabel( 'x positions');

legend( num2str([1:n_link]'-1))
% Numerical
fprintf( 'Numerical stdev at each point : %f %f %f %f %f %f\n', std(x_all,0,1)')
% Predicted
std_pred = [0 sqrt(5/6) sqrt(4/3) sqrt(3/2) sqrt(4/3) sqrt(5/6) ];
fprintf( 'Analytical stdev at each point: %f %f %f %f %f %f\n', std_pred)


subplot(2,1,2);
plot( x_all(:,2), x_all(:,6),'.' );
c = corrcoef(x_all);axis equal
fprintf( 'Corrcoeff between next and previous neighbor: %f . Theoretical: %f\n', c(2,6),0.2);