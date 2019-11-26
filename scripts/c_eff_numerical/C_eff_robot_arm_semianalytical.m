function C_eff = C_eff_semianalytical( sigma, dtheta, n_link, N )
% sigma = std dev. on angle
%  dtheta = fineness of integration for n=4 case.
% n_link = 3 or 4, number of links
% N = order of wrapped Gaussian

if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'N','var') N = 100; end;
if n_link == 4
    if ~exist( 'dtheta','var'); dtheta = 0.0005; end;
    theta = -pi:dtheta:pi;
    tic
    y = wrapped_gaussian(theta,sigma);
    C_eff = 2*pi * sum( y.^4 ./ abs(cos(theta))) * dtheta;
    toc
elseif n_link == 3
    C_eff = 2*pi*2/sqrt(3) * ( wrapped_gaussian(0,sigma)^3 + wrapped_gaussian(2*pi/3,sigma)^3 );
else
    C_eff = NaN;
    fprintf( 'Seminanalytical expression not available yet for n_link = %d\n', n_link );
end
