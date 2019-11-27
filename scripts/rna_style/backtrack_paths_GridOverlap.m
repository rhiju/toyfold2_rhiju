% parameters
NITER = 200000;
sigma = 0.5;
theta0 = 0; %2*pi/10; 
dtheta = 2*pi/10;
dL = 0.1; 
L = 12.0;
n_max= 10;
use_mvksdensity = 1;

%%
% main histogram loop --> forward
hist_forward = {};
for n_link = 1:n_max
    [hist_forward{n_link},pts_forward{n_link},x_range,theta_range] = get_hist_forward(NITER, sigma, theta0, dL, dtheta, L, n_link );
end

%%
% main histogram loop --> reverse
hist_reverse = {};
for n_link = 1:n_max
    [hist_reverse{n_link},pts_reverse{n_link},x_range,theta_range] = get_hist_reverse(NITER, sigma, theta0, dL, dtheta, L, n_link );
end

%%
% plot histograms (x,y), marginalize over theta.
figure(1)
clf
for n_link = 1:10
    subplot(4,5,n_link)
    imagesc( x_range, x_range, sum(hist_forward{n_link},3)' );
    set(gca,'ydir','normal');
    hold on; plot( 0,0,'k.' ); 
    axis image
    axis off
    subplot(4,5,n_link+10);
    imagesc( x_range, x_range, sum(hist_reverse{n_link},3)' );
    set(gca,'ydir','normal');
    hold on; plot( 0,0,'k.' ); 
    axis image
    axis off
end
colormap( 1- copper(100) )
set(gcf, 'PaperPositionMode','auto','color','white');

%%
% C_eff checks
n_link = 10;
for n = 1:(n_link-1)
    C_eff(n) = (2*pi) * sum(hist_forward{n}(:)/NITER .* hist_reverse{n_link-n}(:)/NITER/dtheta/dL^2);
    fprintf( 'C_eff forward %d reverse %d ==> %f\n', n, n_link-n, C_eff(n) );
end


%%
% backtrack to get samples.
figure(2); clf;
n_link = 10;
deltheta = 0.01;
use_mvksdensity = 1;
plot_stuff = 0;

NSAMPLE = 10;
for q = 1:NSAMPLE
    theta = 0;
    x = [0,0];
    xpath = [0,0]';
    tic
    for n = 1:(n_link-2)
        % need to pick theta from a distribution chosen from wrapped
        %  gaussian, but further weighted by reverse histogram.
        possible_theta = [0:deltheta:2*pi];
        p_forward = wrapped_gaussian( (possible_theta-theta)-theta0, sigma );
        
        
        possible_x = x(1) + [cos( possible_theta )];
        possible_y = x(2) + [sin( possible_theta )];
        
        % Issue with straight up interpolation --> many bins go to zero!
        p_reverse = interp3( x_range, x_range, theta_range, hist_reverse{n_link-n},possible_x,possible_y,possible_theta,'linear',0 );
        if use_mvksdensity
            pts = pts_reverse{n_link-n};
            pts(:,3) = mod( pts(:,3), 2*pi );
            s = std( pts )*(4/(3+2)/size(pts,1))^(1/(3+4)); % bandwidth estimator
            xi = [possible_x; possible_y; possible_theta]';
            p_reverse = mvksdensity(pts,xi,'Bandwidth',s)';
        end
        p_reverse = p_reverse/sum(p_reverse)/deltheta;
        
        p = p_forward .* p_reverse;
        p_cumsum = cumsum( p+1e-5 );
        p_cumsum = p_cumsum / p_cumsum(end);
        
        
        theta = interp1( p_cumsum, possible_theta, rand(1) );
        x = x + [cos( theta ),sin(theta)];
        xpath(:,n+1) = x;
        
        if plot_stuff
            subplot(2,1,1);
            plot( possible_theta, [p_forward; p_reverse; p ] );
            subplot(2,1,2)
            imagesc( x_range, x_range, sum(hist_reverse{n_link-n},3)' );
            hold on; plot( possible_x, possible_y, 'k-' );
            axis image
            plot( x(1), x(2), 'kx' );
            pause;
        end
    end
    xpath = [xpath, [-1 0]'];
    xpath = [xpath, [0 0]'];
    toc
    
    
    plot( xpath(1,:), xpath(2,:),'k' );
    hold on; plot( 0,0,'k.' );
    axis image
    axis off
end




