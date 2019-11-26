function x = sample_based_on_neighbors( mu1, std1, mu2, std2 );

mu = (mu1/std1^2 + mu2/std2^2)/ (1/std1^2 + 1/std2^2);
std = sqrt(1/ (1/std1^2 + 1/std2^2));
x = mu + std * randn( 1 );
