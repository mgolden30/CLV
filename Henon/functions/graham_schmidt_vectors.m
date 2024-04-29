function V = graham_schmidt_vectors( x0, n, maxit )
  %{
  PURPOSE:
  Find the singular vectors of J with power iteration, where J is the
  jacobian of f^n.
  %} 
  
  %First generate a trajectory forward in time
  xs = generate_traj( x0, n );

  %Start with an initially orthogonal set of vectors
  V = eye(2);

  for i = 1:maxit
    for j = 1:n
      %fprintf( "%d ", j );
      J = jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J * V;               %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end

    for j = n:-1:1
      %fprintf( "%d ", j );
      J = jacobian( xs(:,j) ); %Jacobian from this point to the next
      V = J.' * V;             %Push vectors forward
      [V, ~] = qr(V);          %Orthogonalize with qr
    end
    %fprintf( "\n" );
  end
end