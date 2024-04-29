function [v_c, v_e] = compute_CLVs( x0, n, maxit )
  %{
  PURPOSE:
  Compute the CLVs

  v_c - contracting
  v_e - expanding
  %}
  
  %% Contracting Direction

  %First generate a trajectory forward in time
  xs = generate_traj( x0, n );
  V = eye(2); %Start with an initially orthogonal set of vectors
  for i = 1:maxit
    for j = 1:n
      J = jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J * V;               %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end

    for j = n:-1:1
      J = jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J.' * V;             %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end
  end
  v_c = V(:,2);

  %% Now do expanding
  xs = generate_backwards_traj( x0, n );
  V = eye(2); %Start with an initially orthogonal set of vectors
  for i = 1:maxit
    for j = 1:n
      J = backwards_jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J * V;               %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end

    for j = n:-1:1
      J = backwards_jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J.' * V;             %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end
  end
  v_e = V(:,2);
end