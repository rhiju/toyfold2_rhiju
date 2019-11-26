function hist_forward = update_histogram( x, theta, hist_forward, L, dL, dtheta );
bin1 = round( (x(1) + L)/dL ) + 1;
bin2 = round( (x(2) + L)/dL ) + 1;
Ntheta = size( hist_forward, 3);
Nxy = size( hist_forward, 1);
bin3 = mod( round( theta/dtheta ) - 1, Ntheta ) + 1;
if (bin1 > 0 & bin1 <= Nxy) & (bin2 > 0 & bin2 <= Nxy )
    hist_forward(bin1,bin2,bin3) = hist_forward(bin1,bin2,bin3)+1;
end
