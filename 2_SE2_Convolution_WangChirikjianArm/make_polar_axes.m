hold on
p = [0:0.01:2*pi];
plot3(cos(p),sin(p),    0*p+1,'color',[0.5 0.5 0.5] );
plot3(2*cos(p),2*sin(p),0*p+1,'color',[0.5 0.5 0.5] );
plot3(3*cos(p),3*sin(p),0*p+1,'color',[0.5 0.5 0.5] );
plot3(4*cos(p),4*sin(p), 0*p+1,'color',[0.5 0.5 0.5] );
for i = 1:12;
    plot3( 4*cos(2*pi*i/12) * [-1 1], 4*sin(2*pi*i/12)*[-1 1], [1 1], 'color',[0.5 0.5 0.5]);
end
axis( 4*[-1 1 -1 1]);
axis equal
