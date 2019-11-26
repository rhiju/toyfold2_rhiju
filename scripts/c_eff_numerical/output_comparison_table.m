function output_comparison_table( sigma_test, titles, sigma_val, Ceff_val );
% output_comparison_table( sigma_test, titles, sigma_val, Ceff_val );

for n = 1:length(titles);    
    if ( n > 1 )fprintf( '\t' );  end;
    fprintf( 'Rhiju %s',titles{n} );
end
fprintf( '\n' );
for k = 1:length(sigma_test );
    for n = 1:length(titles);
        if ( n > 1 )fprintf( '\t' );  end;
        % average over any values that have the same sigma...
        sigma_unique = unique( sigma_val{n} );
        Ceff_unique = [];
        for i = 1:length(sigma_unique); Ceff_unique(i) = mean( Ceff_val{n}( sigma_val{n} == sigma_unique(i) ) ); end;
        val = exp( interp1( sigma_unique,log( Ceff_unique ),sigma_test(k) ) );
        fprintf( '%f',val );        
    end
    fprintf( '\n' );
end
fprintf( '\n' );




