N = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open N-linked chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros( N, N );
A(1,1) = 1;
A(1,2) = -1;
A(2,2) = -1;
for i = 2:N-1
    A(i,i) = 2;
    A(i,i-1) = -1;
    A(i,i+1) = -1;
    A(i-1,i) = -1;
    A(i+1,i) = -1;
end
A(N,N) = 1;
A(N,N-1) = -1;
A(N-1,N) = -1;
[V,D] = eig(A);
D = diag(D);
D'
sort( 4*sin(pi*[0:1/2:(N-1)/2]/N).^2 )
fprintf( '\nProduct of non-zero eigenvalues OPEN  CHAIN: %f\n', prod(D(2:end)) )

clf
subplot(2,2,1);
imagesc( V' );
title( 'open chain' ); ylabel( 'Mode' ); xlabel( 'Bead position');

subplot(2,2,2);
title( 'open chain -- analytical vs. numerical' )
X = [0.5:0.1:N+0.5];
colors = hsv(N);
for q = 1:N
    hold on
    plot( V(:,q)/sign(V(1,q)), 'o','markeredgecolor',colors(q,:),'markerfacecolor',colors(q,:) );
    v = cos( 2*pi*(X-0.5)*(q-1)/2/N )/    sqrt(sum(cos( ([1:N]-0.5)*2*pi*(q-1)/2/N ).^2));
    v = v/sign(v(X==1.0));
    plot( X,v, '-','color',colors(q,:) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After chain closure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 2 * eye(N) + circshift(eye(N),[0,1]) + circshift(eye(N),[0,-1]);
[V,D] = eig(A);
D = diag(D);
D'
sort( 4*sin(pi*[1:N]/N).^2 )
fprintf( 'Product of non-zero eigenvalues CLOSED CHAIN: %f\n', prod(D(2:end)) )

subplot(2,2,3);
imagesc( V' );
title( 'closed chain' ); ylabel( 'Mode' ); xlabel( 'Bead position');

X = [1:0.1:6];
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(2,2,4);
hold on;
X = [0:0.1:N];
for q = [N,1:N-1]
    hold on
    plot( V(:,q)/sign(V(1,q)), 'o','markeredgecolor',colors(q,:),'markerfacecolor',colors(q,:) );
    phase_shift = 0;
    v = cos( 2*pi*X*q/N + phase_shift)/ sqrt(sum(cos( ([1:N])*2*pi*q/N +phase_shift).^2));
    v = v/sign(v(X==1.0));
    plot( X,v, '-','color',colors(q,:) );
end
title( 'closed chain -- analytical vs. numerical [up to phase shift]' ); ylabel( 'Mode' ); xlabel( 'Bead position');


