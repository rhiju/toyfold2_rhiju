function C_eff = C_eff_compute_MM2( sigma );

% theta = value to evaluate function at 
% Sigma = std dev. of wrapped gaussian
% N = total number of wraps (ideally should be infinity)
N = 500;
theta_all = 2*pi*[-N:N];

Z = sum( exp( -theta_all.^2/2/sigma^2 ) );
y = 1/sqrt(2*pi)/sigma * Z;
soft_factor = 1 - sum(theta_all.^2/sigma^2 .* exp( -theta_all.^2/2/sigma^2 ) )/Z;
sigma_eff = sigma/sqrt(soft_factor);
if soft_factor < 0; C_eff = NaN; return; end;
C_eff = (2*pi)^(3/2) * y^4 * (sigma_eff/2)* erf( sqrt(2) * pi/sigma_eff);