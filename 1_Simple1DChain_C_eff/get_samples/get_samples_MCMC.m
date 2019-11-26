NITER = 100000;
n_link = 6;
x = [0, 0.1*randn( 1, n_link-1 )]; % deviation of each link position from its ideal value.
sep_ideal = [-1,1,-1,1,-1,1];
x = x + [0,1,0,1,0,1];
% metropolis chain monte carlo
current_score = chain_score(x, sep_ideal);
step_size = 0.2;
x_all = zeros(NITER,n_link);
score_all = zeros(1,NITER);
accepts = 0;
for i = 1:NITER    
    x_new = x;
    move_idx = 1 + randi(5);
    x_new( move_idx ) =  x_new( move_idx ) + step_size * randn(1);
    
    new_score = chain_score(x_new, sep_ideal);
    if new_score < current_score || ...
       exp(  current_score - new_score  ) > rand(1)   
       accepts = accepts + 1;
       x = x_new;
       current_score = new_score;
       x_all(i,:) = x;
       score_all(i) = current_score;
    end
end
plot(x_all);
accepts/NITER

legend( num2str([1:n_link]'-1))
% Numerical
std(x_all,0,1)
% Predicted
[0 sqrt(5/6) sqrt(4/3) sqrt(3/2) sqrt(4/3) sqrt(5/6) ]