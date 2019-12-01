function pts = get_pts_forward( NITER, sigma, theta0, n_forward, trans,rot,roadblock )
% pts = get_pts_forward( NITER, sigma, theta0, n_forward )
%
% Get sampled x,y,theta for n-link polymer
%
% INPUTS:
%  NITER  = number of independent samples
%  sigma  = width of wrapped Gaussian for angle (radians)
%  theta0 = preferred angle at each link (radians)
%  n_forward = number of polymer links
%
% OPTIONAL INPUTS:
%  trans    = (x,y) start translation  [default 0,0]
%  rot      = start rotation (radians) [default 0]
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
    while ( n <= n_forward )
        theta_next = theta + theta0 + sigma*randn(1);
        x_next = x + [cos( theta_next) sin(theta_next) ];
        if exist( 'roadblock', 'var' ) && ~isempty( roadblock )
            if norm( x_next - roadblock(1:2) ) < roadblock(3); continue; end;
        end            
        n = n+1;
        theta = theta_next;
        x = x_next;
        xpath(:,n) = x;
    end
    pts = [pts; x,theta];
    if (NITER <= 100 ); plot( xpath(1,:),xpath(2,:), '-','linew',0.5,'color',[0.5,0.5,0.5] ); hold on; end;
end

