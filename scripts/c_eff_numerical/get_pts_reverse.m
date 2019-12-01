function pts = get_pts_reverse( NITER, sigma, theta0, n_reverse, trans,rot,roadblock )
% pts = get_pts_reverse( NITER, sigma, theta0, n_reverse )
%
% Get sampled x,y,theta for n-link polymer
%
% INPUTS:
%  NITER  = number of independent samples
%  sigma  = width of wrapped Gaussian for angle (radians)
%  theta0 = preferred angle at each link (radians)
%  n_reverse = number of polymer links
%  roadblock = [ctr_x,ctr_y,radius] of circular steric exclusion zone.
%
% OPTIONAL INPUTS:
%  trans    = (x,y) target translation  [default 0,0]
%  rot      = target rotation (radians) [default 0]
%  roadblock = [ctr_x,ctr_y,radius] of circular steric exclusion zone.
%
% OUTPUTS:
%  pts = NITER samples of (x, y, theta) 
%
% (C) R. Das, Stanford University 2019
if ~exist( 'trans', 'var') | length(trans) < 2; trans = [0,0]; end;
if ~exist( 'rot', 'var') rot = 0; end;

pts = [];
for i = 1:NITER
    x = trans;
    theta = rot;
    n = 1;
    xpath(:,n) = x';
    while ( n <= n_reverse )
        x_prev = x - [cos( theta) sin(theta) ];
        theta_prev = theta - theta0 - sigma*randn(1);
        if exist( 'roadblock', 'var' )  && ~isempty( roadblock )
            % need to predict one step ahead...
            x_prev2 = x_prev - [cos(theta_prev) sin(theta_prev)]; 
            if norm( x_prev2 - roadblock(1:2) ) < roadblock(3); continue; end;
        end            
        n = n+1;
        theta = theta_prev;
        x = x_prev;
        xpath(:,n) = x;
    end
    pts = [pts; x,theta];
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
end