function s = chain_score( x, sep_ideal );

s = 0.5 * sum( (x - circshift(x,1) - sep_ideal).^2 );
