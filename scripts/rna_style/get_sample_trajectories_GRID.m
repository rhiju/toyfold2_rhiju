function all_xpaths = get_sample_trajectories_GRID( NSAMPLE, all_n_link, all_trans, trans_vec, rot, theta0, sigma, pts_reverse,  use_mvksdensity, plot_steps );

clf;
count = 0;
for i = 1:length(all_n_link)
    n_link = all_n_link(i);
    for j = 1:length(all_trans);
        count = count+1;
        trans = trans_vec*all_trans(j);
        subplot( length(all_n_link), length(all_trans), count );
        all_xpath = get_sample_trajectories( NSAMPLE, n_link, theta0, sigma, pts_reverse, trans,rot, use_mvksdensity, plot_steps);
        all_xpaths{i,j} = all_xpath;
        drawnow;
        set(gca,'ylim',[-9 9]);
        if ( j == 1 ); set(title( ['n_{link} = ',num2str(n_link)] ),'interpreter','tex'); end;
    end
end




