function log_y = log_wrapped_gaussian( theta, sigma, N);
%l og_y = log_wrapped_gaussian( theta, sigma, N);
%
% Slightly different funciton than wrapped_gaussian -- to handle when
%  gaussian numerically evaluates to zero.
%
% theta = value to evaluate function at -- can be multiple values.
% Sigma = std dev. of wrapped gaussian
% N = total number of wraps (ideally should be infinity)

if ~exist( 'N', 'var'); N = 100; end;
theta_all = repmat(theta, 2*N+1, 1) + repmat( 2*pi*[-N:N]', 1, length(theta) );

contribs = -log(sqrt(2*pi))-log(sigma)-theta_all.^2/2/sigma^2;
log_y = max( contribs ) + log( sum( exp( contribs - max( contribs ) ) ) );
    
    


