function f = lyapunov_objective(L, zs, n, T)
  %{
  PURPOSE:
  Find the matrix logarithm of the monodromy matrix.
  %}

  L = reshape( L, [3,3] );
  dt = T/n;

  phi = eye(3);
  for i = 1:n
    %Use Crank-Nicolson to solve for the next timestep
    J  = jacobian( zs(:,i)   );
    Jp = jacobian( zs(:,i+1) );

    %phi = phi + dt/2*( J*phi - phi*L );
    %phi = sylvester( eye(3)/2 - dt/2*Jp, eye(3)/2 + dt/2*L, phi );
    phi = sylvester( eye(3)/2 - dt*Jp, eye(3) + dt*L, phi );
  end

  f = phi - eye(3);
  f = reshape(f, [9,1] );

  v = rossler(zs(:,1));
  f = [f; L*v];
end