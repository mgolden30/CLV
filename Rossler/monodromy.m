%{
Converge monodromy matrices for periodic orbits
%}

clear;


addpath("functions");
load("POs/3.mat");

hold on
zs = rk4( z(1:3), z(4)/n, n);
scatter3( zs(1,:), zs(2,:), zs(3,:) );

phi = estimate_monodromy( zs, n, z(4) )

e = eigs(phi)

e(3)
return


%L = zeros(3,3);
L = reshape(L, [9,1]);

h = 1e-3;
maxit = 32;
for i =1:maxit
  f = lyapunov_objective(L, zs, n, z(4) );
  
  J = zeros(12,9);
  for j = 1:9
    L2 = L;
    L2(j) = L2(j) + h;
    J(:,j) = ( lyapunov_objective(L2, zs, n, z(4) ) - f)/h;
  end

  tol = 1;
  [dL, ~] = lsqr(J, f, tol);
  L = L - dL;

  %At each step, enforce orthogonality to the starting velocity
  %v = rossler(zs(:,1));
  %v = v/norm(v);
  %P = eye(3) - v*v';
  %L = reshape(L,[3,3]);
  %L = L*P;
  %L = reshape(L,[9,1]);

  norm(f)
end

%%
L = reshape(L, [3,3])
trace(L)
eigs(L)

v = rossler(zs(:,1));

(L*v)./v
