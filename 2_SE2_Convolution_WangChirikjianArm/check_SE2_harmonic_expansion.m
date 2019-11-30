figure(2)
clear i

% Parameters of link arm distribution
theta0 = 30 * pi/180;  % mean turn angle
sigma = sqrt(0.1);     % standard deviation in turn angle (wrapped gaussian)

% How many orders to keep in Fourier transform 
delp = 0.2;
p_max = 300;
p_val = [0:delp:p_max];
max_order = 10;
m_val = -max_order:1:max_order;
n_val = -max_order:1:max_order;

delr = 0.05;
r_val = [0:delr:8];
del_phi = 0.1;
phi_val = [-pi:del_phi:pi+0.05];

x_val   = [-4:0.1:4];
y_val   = [-4:0.1:4];


% show 1 harmonic basis function:
%n = 1;
%m = 1;
%p = 1;
% u = i^(n-m) * exp( - i * (n*thetag + (m-n)*phig )) .* besselj( n-m, p * rg );        
%clf
%imagesc( real( u(:,:,1) ) )

tic

rho_hat_final = zeros(length(m_val),length(n_val), length(p_val) );
rho_hat = zeros(length(m_val),length(n_val) );
fprintf( '\n\nComputing transform of 3 arms...\n' );
for p_idx = 1:length(p_val)
    p = p_val(p_idx);
    % could accelerate following through ndgrid.
    for m_idx = 1:length(m_val)
        for n_idx = 1:length(n_val)
            m = m_val(m_idx);
            n = n_val(n_idx);
            rho_hat(m_idx,n_idx) = i^(n-m) * besselj( m-n, p ) * exp( -n*n*sigma*sigma ) * exp( i * n * theta0 );
        end
    end
    % the convolution becomes a straight-up matrix multiplication
    rho_hat_final(:,:,p_idx) = rho_hat * rho_hat * rho_hat;
end
toc

%m_idx = find( m_val == 0 );
%imagesc( squeeze(real(rho_hat_final(m_idx,:,:) ) ) );
%return;

% now transform back to real-space
% do average over theta (so just plotting location of end of robot arm).
% need to use rho_final(m,0,p)  to get c_m( r )
tic
fprintf( '\nTransform back to polar coordinates, step 1 (integral over p)...\n' );
n_idx = find( n_val == 0 );
c = zeros( length(m_val), length(r_val) );
for r_idx = 1:length(r_val )
    r = r_val(r_idx);
    for m_idx = 1:length(m_val)
        m = m_val(m_idx);
        % which one of the following is it?
        % what I got...
        %c(m_idx,r_idx) = sum( squeeze(rho_hat_final(m_idx,n_idx,:))' .* besselj( -m, p_val*r ) .* p_val * delp);
        % in wang paper
        c(m_idx,r_idx) = sum( squeeze(rho_hat_final(n_idx,m_idx,:))' .* besselj( -m, p_val*r ) .* p_val * delp);
    end
end
toc

% and now work out rho( r, phi );
tic
fprintf( 'Transform back to polar coordinates, step 2 (integral over harmonic order m)...\n' );
rho_final = zeros(length(r_val),length(phi_val) );
for r_idx = 1:length(r_val )
    r = r_val(r_idx);
    for phi_idx = 1:length( phi_val )
        phi = phi_val(phi_idx);
        rho_final(r_idx,phi_idx) = sum( i.^-m_val .* exp( -i * m_val * phi ) .* c(:,r_idx)' * del_phi);
    end
end
toc

clf
%imagesc( r_val,phi_val,real(rho_final') );
%xlabel( 'r'); ylabel( 'phi')

tic
%fprintf( '\nDoing polar-to-cartesian interpolation...\n' );
%for x_idx = 1:length(x_val )
%     for y_idx = 1:length(y_val )
%         x = x_val(x_idx );
%         y = y_val(y_idx );
%         r = sqrt( x*x + y*y);
%         phi = atan2( y, x );
%         rho_xy( x_idx, y_idx ) = interpn( r_grid, phi_grid, rho_final, r, phi);
%     end
% end


clf
[r_grid, phi_grid ] = ndgrid( r_val, phi_val );
[x_grid, y_grid ] = pol2cart( phi_grid, r_grid );
surf( x_grid, y_grid, real(rho_final) )
shading interp
view(2)
axis equal
toc

%subplot(2,1,2)
%imagesc( x_val,y_val,real(rho_xy') );
xlabel( 'x'); ylabel( 'y')
set(gca,'ydir','normal');
colormap( 1-copper(100));
make_polar_axes
fprintf( 'Maximum value: %f    Minimum value: %f\n',  max(max(real(rho_xy))), min(min(real(rho_xy))) )

