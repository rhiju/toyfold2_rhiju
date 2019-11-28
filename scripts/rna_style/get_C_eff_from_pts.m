function C_eff = get_C_eff_from_pts( pts_f, pts_r);
% C_eff = get_C_eff_from_pts( pts_f, pts_r);

pts_r(:,3) = mod( pts_r(:,3), 2*pi );
pts_f(:,3) = mod( pts_f(:,3), 2*pi );

% s = get_kde_bandwidth( pts_f );
% p = mvksdensity(pts_f,pts_r,'Bandwidth',s)';

% external library that fails (and is slower too)
%p = akde(pts_f,pts_r); %,'Bandwidth',s)';

% try to slice out theta ranges, and do (x,y) KDE's separately
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