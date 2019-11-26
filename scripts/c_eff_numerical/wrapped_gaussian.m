function y = wrapped_gaussian( theta, sigma, N);
% theta = value to evaluate function at -- can be multiple values.
% Sigma = std dev. of wrapped gaussian
% N = total number of wraps (ideally should be infinity)
if ~exist( 'N', 'var'); N = 100; end;
theta_all = repmat(theta, 2*N+1, 1) + repmat( 2*pi*[-N:N]', 1, length(theta) );
y = 1/sqrt(2*pi)/sigma * sum(  exp( -theta_all.^2/2/sigma^2 ) );

