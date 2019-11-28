function all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans, rot, use_mvksdensity, plot_steps);
% all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_forward, pts_reverse, use_mvksdensity, plot_steps);
%
% Sample trajectories based on forward steps, weighted by reverse
%   probability density (using KDE by default)
%
%
% (C) R. Das, Stanford University

for n = 1:(n_link-1)
    pts_reverse{n}(:,1:2) = pts_reverse{n}(:,1:2) * [cos(rot), sin(rot); -sin(rot) cos(rot)] + trans;
    pts_reverse{n}(:,3) = pts_reverse{n}(:,3) + rot;
end

if plot_steps | ~use_mvksdensity
    tic
    dtheta = 2*pi/10;dL = 1; L = 12.0;
    for n = 1:(n_link-1)
        [hist_reverse{n},x_range,theta_range] = get_histogram(pts_reverse{n}, L, dL, dtheta);
    end
    toc
end

deltheta = 0.01;
all_xpath = zeros(2,n_link+1,NSAMPLE);
for q = 1:NSAMPLE
    fprintf( 'Evaluating trajectory %d of %d...\n', q, NSAMPLE );
    theta = 0;
    x = [0,0];
    xpath = [0,0]';
    tic
    for n = 1:(n_link-1)
        % need to pick theta from a distribution chosen from wrapped
        %  gaussian, but further weighted by reverse histogram.
        possible_theta = [0:deltheta:2*pi];

        % One-dimensional probability density (forward):
        p_forward = wrapped_gaussian( (possible_theta-theta)-theta0, sigma );
        
        possible_x = x(1) + [cos( possible_theta )];
        possible_y = x(2) + [sin( possible_theta )];
        
        % One-dimensional probability density (reverse)
        if use_mvksdensity
            pts = pts_reverse{n_link-n};
            pts(:,3) = mod( pts(:,3), 2*pi );
            xi = [possible_x; possible_y; possible_theta]';

            % MATLAB's multivariate ksdensity toolbox. 
            s = get_kde_bandwidth( pts );
            p_reverse = mvksdensity(pts,xi,'Bandwidth',s)';
            
            % Tried this other toolbox (circ_ksdensityn) but it appears to stall?
            %p_reverse = circ_ksdensityn(pts(1:10000,:),xi,[NaN NaN; NaN NaN; 0 2*pi])';
        else
            % Issue with straight up interpolation --> many bins go to zero!
            p_reverse = interp3( x_range, x_range, theta_range, hist_reverse{n_link-n},possible_x,possible_y,possible_theta,'linear',0 );
        end
        if ( sum(p_reverse) == 0.0 );  xpath(1:2,[n:n_link] ) = NaN; break; end;
        p_reverse = p_reverse/sum(p_reverse)/deltheta;
        
        p = p_forward .* p_reverse;
        p_cumsum = cumsum( p+1e-10 );
        p_cumsum = p_cumsum / p_cumsum(end);

        if ( n < n_link-1 )
            theta = interp1( p_cumsum, possible_theta, rand(1) );
            x = x + [cos( theta ),sin(theta)];
        else
            % penultimate point is specified at a special spot -- 
            %  one length 'back' from endpoint.
            x = [-1, 0] * [cos(rot), sin(rot); -sin(rot) cos(rot)] + trans;
        end
        
        xpath(:,n+1) = x;
        if plot_steps
            subplot(2,1,1);
            plot( possible_theta, [p_forward; p_reverse; p ],'linew',2 );
            set(legend( 'p_forward','p_reverse','p' ),'interp','none');
            subplot(2,1,2)
            imagesc( x_range, x_range, sum(hist_reverse{n_link-n},3)' );
            hold on; plot( possible_x, possible_y, 'k-' );
            axis image
            plot( x(1), x(2), 'kx' );
            plot( xpath(1,:),xpath(2,:),'k' );
            colormap( 1 - copper(100) );
            legend( 'Possible next positions','Choice' );
            set(gcf, 'PaperPositionMode','auto','color','white');
            set(gca,'ydir','normal');
            pause;
        end
    end
    
    % last points are at fixed locations:
    xpath = [xpath, trans'];
    toc    

    all_xpath(:, :, q) = xpath;
end


cla;
for n = 1:NSAMPLE
    plot( all_xpath(1,:,n), all_xpath(2,:,n),'linewidth',2 ); hold on
end
plot( 0,0,'ko','markerfacecolor','k','markersize',5 ); 
plot( trans(1),trans(2),'kx','markersize',5 ); 
axis image; axis off



% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = get_transform_matrix( trans, rot )
% 
% M = [cos( rot ) -sin( rot ) trans(1); sin(rot), cos(rot) trans(2); 0 0 1]';
% M
% 
