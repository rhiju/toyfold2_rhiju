function C_eff = get_C_eff_from_pts( pts_f, pts_r);
% C_eff = get_C_eff_from_pts( pts_f, pts_r);
%
% INPUTS
%  pts_f = sampled (x,y,theta) for forward chain segment
%  pts_r = sampled (x,y,theta) for reverse chain segment
%
% OUTPUT
%  C_eff = Effective molarity for chain closure between forward samples
%            and reverse samples, calculated based on overlap of
%            distributions (over translations & rotations), and using
%            MATLAB's KDE estimator.
%
% (C) R. Das, Stanford University, 2019

% need to collapse angles into 0 to 2*pi range.
pts_r(:,3) = mod( pts_r(:,3), 2*pi );
pts_f(:,3) = mod( pts_f(:,3), 2*pi );

s = get_kde_bandwidth( pts_f );
p = mvksdensity(pts_f,pts_r,'Bandwidth',s)';

% external library that fails (and is slower too)
%p = akde(pts_f,pts_r); %,'Bandwidth',s)';

% also tried a circular KDE library avaliable on Mathwork File Exchange
%   but kept getting garbage.

% try to slice out theta ranges, and do (x,y) KDE's separately
%  fails to get good estimates -- turns out to be really powerful to have
%  KDE in theta direction.
% dtheta = pi/10;
% theta_range = [0:dtheta:2*pi];
% p = [];
% for i = 1:length(theta_range)-1
%     idx_f = find( pts_f(:,3)>= theta_range(i) & pts_f(:,3) < theta_range(i+1) );
%     idx_r = find( pts_r(:,3)>= theta_range(i) & pts_r(:,3) < theta_range(i+1) );
%     if ( length( idx_r ) == 0 | length( idx_f) == 0 ) continue; end;
%     s = get_kde_bandwidth( pts_f(idx_f,1:2) );
%     p = [p, mvksdensity(pts_f(idx_f,1:2),pts_r(idx_r,1:2),'Bandwidth',s)'/dtheta];
%     %s = get_kde_bandwidth( pts_f );
%     %p = [p, mvksdensity(pts_f,pts_r(idx_r,:),'Bandwidth',s)'];
% end

C_eff = (2*pi) * mean(p);
