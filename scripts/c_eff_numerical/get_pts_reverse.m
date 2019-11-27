function pts = get_pts_reverse( NITER, sigma, theta0, n_reverse )
% pts = get_pts_reverse( NITER, sigma, theta0, n_reverse )
%
% Get sampled x,y,theta for n-link polymer
%
% INPUTS:
%  NITER  = number of independent samples
%  sigma  = width of wrapped Gaussian for angle (radians)
%  theta0 = preferred angle at each link (radians)
%  n_reverse = number of polymer links
%
% OUTPUTS:
%  pts = NITER samples of (x, y, theta) 
%
% (C) R. Das, Stanford University 2019

pts = [];
for i = 1:NITER
    theta = 0;
    x = [0,0];
    xpath = [0,0]';
    for n = 1:n_reverse
        x = x - [cos( theta) sin(theta) ];
        theta = theta - theta0;
        theta = theta - sigma*randn(1);
        xpath(:,n+1) = x;
    end
    pts = [pts; x,theta];
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
end