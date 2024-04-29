%{
Converge periodic orbits variationally
%}

addpath("functions/");

x0 = [1;2;3];

n = 1024*3;
dt = 0.025;
xs = rk4( x0, dt, n );

tr = 256; %transient
ms = 10;
scatter3( xs(1,tr:end), xs(2,tr:end), xs(3,tr:end), ms, 'filled' );
%return;

%Let's find a cycle or two
l = 3; %cycle length
T = 5.88*l;
n = 256*l;
zs = [xs(:,end-(n)+1:end)];
zs = reshape(zs, [3*n,1]);
zs = [zs; T];

%%
h = 1e-4;
maxit = 256;
for i =1:maxit
  
  f = variational_objective( zs );
  fprintf("%d: |f| = %e\n", i, norm(f) );

  J = @(v) (variational_objective( zs + h*v ) - f)/h;
  
  inner = 64;
  outer = 1;
  tol   = 1e-3;
  hook  = 0.1;

  [dz, ~] = gmres(J, f, inner, tol, outer);
  zs = zs - hook*dz;
end

traj = reshape( zs(1:end-1), [3,n] );
scatter3( traj(1,:), traj(2,:), traj(3,:) );

return;

%% Save POs
save("variational_POs/4.mat", "zs", "n");

function f = periodic_objective(z,n)
  x0 = z(1:3);
  T  = z(4);

  %Shoot forward
  dt = T/n;
  xs = rk4( x0, dt, n );
  
  f = xs(:,end) - x0;
end


