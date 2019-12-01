d0 = load( 'save_roadblockNONE_50link.mat', 'all_C_eff' );
d1 = load( 'save_roadblock10_50link.mat', 'all_C_eff' );
n_max = size( d0.all_C_eff, 2 );
clf
for n = 6:n_max
    plot( n, d0.all_C_eff(3:(n-3),n),'bo','markerfacecolor','b' ); hold on
    plot( n, d1.all_C_eff(3:(n-3),n),'ro','markerfacecolor','r' );
end
xlabel( 'number of links' )
ylabel( 'C_{eff}' );
legend( 'without roadblock','with roadblock' );
set(gcf, 'PaperPositionMode','auto','color','white');
title( 'Grid-based estimation of C_{eff}' );



