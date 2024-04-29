%{
Converge monodromy matrices for periodic orbits with variational methods
%}

clear;

addpath("functions");

%Load an orbit of your choice
load("variational_POs/1.mat");

T = zs(end);
traj = reshape( zs(1:end-1), [3,n] );

[phi, monodromy]  = estimate_monodromy( traj, n, T );

L = logm( monodromy )/T;
eigs(L)

for i = 1:n
  phi(:,:,i) = eye(3);
end

plot(squeeze(abs(phi(1,1,:))));

s = pack_state(phi, L);

%%
h = 1e-4;
maxit = 128;

for i =1:maxit
  f = variational_objective_monodromy( s, traj, T );
  J = @(v) ( variational_objective_monodromy( s+h*v, traj, T ) - f)/h;
  
  fprintf("%d: |f| = %e\n", i, norm(f) );
  inner = 64;
  outer = 1;
  tol   = 1e-9;
  hook  = 0.1;

  [ds, ~] = gmres(J, f, inner, tol, outer);
  s = s - hook*ds;
  s = optimize_L( s, traj, T ); 
end

[phi,L] = unpack_state(s);

eigs(L)
