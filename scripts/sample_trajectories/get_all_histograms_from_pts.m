function [hist_forward, hist_reverse,x_range,theta_range] = get_all_histograms_from_pts( pts_forward, pts_reverse, dL, dtheta, L );
clf;    
nrows = ceil(length(pts_forward)/5);
for n = 1:length(pts_forward);
    subplot(2*nrows,5,n)
    [hist_forward{n},x_range,theta_range] = get_histogram(pts_forward{n}, L, dL, dtheta);
    imagesc( x_range, x_range, sum(hist_forward{n},3)' );
    title( sprintf('Forward: %d links',n) );
    set(gca,'ydir','normal');
    hold on; plot( 0,0,'k.' ); 
    axis image
    axis off

    subplot(2*nrows,5,n+nrows*5);
    [hist_reverse{n},x_range,theta_range] = get_histogram(pts_reverse{n}, L, dL, dtheta);
    imagesc( x_range, x_range, sum(hist_reverse{n},3)' );
    title( sprintf('Reverse: %d links',n) );
    set(gca,'ydir','normal');
    hold on; plot( 0,0,'k.' ); 
    axis image
    axis off
end
colormap( 1- copper(100) )
set(gcf, 'PaperPositionMode','auto','color','white');
