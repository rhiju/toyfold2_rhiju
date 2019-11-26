function [C_eff,rho_final] = C_eff_SE2_harmonic_expansion( sigma, p_max, max_order, n_link, theta0 );
% [C_eff,rho_final] = C_eff_SE2_harmonic_expansion( sigma, p_max, max_order, n_link, theta0 );

clear i

% Parameters of link arm distribution
if ~exist( 'n_link','var') n_link = 4; end;
if ~exist( 'p_max','var') p_max = 300; end;
if ~exist( 'max_order','var') max_order = 10; end;
if ~exist( 'theta0','var') theta0 = 2*pi/n_link; end;

% How many orders to keep in Fourier transform 
delp = 0.2;
p_val = [0:delp:p_max];
m_val = -max_order:1:max_order;
n_val = -max_order:1:max_order;

delr = 0.05;
r_val = [0:delr:8];
del_phi = 0.1;
phi_val = [-pi:del_phi:pi+0.05];

x_val   = [-4:0.1:4];
y_val   = [-4:0.1:4];

tic
rho_hat_final = zeros(length(m_val),length(n_val), length(p_val) );
rho_hat = zeros(length(m_val),length(n_val) );
fprintf( '\n\nComputing transform of %d arms...\n', n_link );
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
    rho_hat_final(:,:,p_idx) = rho_hat;
    for k = 1: (n_link-1)
        rho_hat_final(:,:,p_idx) = rho_hat_final(:,:,p_idx) * rho_hat;
    end
end
toc

% to get value of probability density (d^3 p/dxdydtheta), need to evaluate
% rho at r = 0, phi = 0, theta = 0, which ends up being a neat trace over
% rho_hat, as follows
rho_hat_final_sum = 0.0;
for p_idx = 1:length(p_val)
    p = p_val(p_idx);
    for m_idx = 1:length(m_val)
        m = m_val(m_idx);
        rho_hat_final_sum = rho_hat_final_sum  + rho_hat_final(m_idx,m_idx,p_idx)*p*delp;
    end
end

% the 2 pi converts to effective concentration 
% (normalize to standard state with (1/2 pi) rad^-1).

% Added a couple more 1/(2*pi)^2 based on what I think might be missing in
% Chirijkian-Wang equations. Those would also normally appear in a 2D
% Fourier transform, and make the numberscome out nice.
C_eff = 2*pi * rho_hat_final_sum/(2*pi)^2;


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
rho_final_sum = 0;
for r_idx = 1:length(r_val )
    r = r_val(r_idx);
    for phi_idx = 1:length( phi_val )
        phi = phi_val(phi_idx);
        rho_final(r_idx,phi_idx) = sum( i.^-m_val .* exp( -i * m_val * phi ) .* c(:,r_idx)' * del_phi);
        rho_final_sum = rho_final_sum + rho_final(r_idx,phi_idx)*delr*r*del_phi;
    end
end

rho_final_sum
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
max(max(rho_final) )
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
%fprintf( 'Maximum value: %f    Minimum value: %f\n',  max(max(real(rho_xy))), min(min(real(rho_xy))) )

