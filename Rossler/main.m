x0 = [1;2;3];
%{
x0 = [
   2.6294
   -5.3157
    0.0481    
];
%}
n = 1024*2;
dt = 0.025;
xs = rk4( x0, dt, n );

tr = 256; %transient
ms = 10;
scatter3( xs(1,tr:end), xs(2,tr:end), xs(3,tr:end), ms, 'filled' );
%return;


%Let's find a cycle or two
l = 5; %cycle length
T = 5.88*l;
z = [xs(:,end); T];
n = 256*l;

h = 1e-4;
maxit = 512;
for i =1:maxit
  f = periodic_objective(z, n);

  J = zeros(3,4);
  for j = 1:4
    z2 = z;
    z2(j) = z2(j) + h;
    J(:,j) = (periodic_objective(z2,n) - f)/h;
  end

  [dz, ~] = lsqr(J, f);
  z = z - 0.1*dz;
  norm(f)
end


z'

zs = rk4( z(1:3), z(4)/n, n);
scatter3( zs(1,:), zs(2,:), zs(3,:) );


%% Save POs
save("POs/3.mat", "z", "n");

function f = periodic_objective(z,n)
  x0 = z(1:3);
  T  = z(4);

  %Shoot forward
  dt = T/n;
  xs = rk4( x0, dt, n );
  
  f = xs(:,end) - x0;
end


