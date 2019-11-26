NITER = 1000000;
C_eff_sum = 0.0;
n_link = 2;
for i = 1:NITER
    x = randn( 1, n_link-1 ); % deviation of each link position from its ideal value.
    y = sum( x );
    C_eff_sum = C_eff_sum + (1/sqrt(2*pi)) * exp( -y^2/2);
end

C_eff_numerical = C_eff_sum/NITER;
C_eff_exact = (1/sqrt(2*pi*n_link));
fprintf( 'Number of links: %d ==> Numerical: %f  Exact: %f \n', n_link, C_eff_numerical, C_eff_exact )


%%
% Harmonic mode analysis -- not quite there.

for N = 1 :10
    A = zeros(N);
    for i = 1:N;
        A(i,i) = 2;
        A(i,mod(i+1 -1,N)+1) = -1;
        A(i,mod(i-1 -1,N)+1) = -1;
    end;
    [V,D ] = eig(A);
    dvals = diag(D)
    plot(V)
    sqrt(prod(dvals(2:end)))
end




