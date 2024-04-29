function f = variational_objective_monodromy( s, traj, T )
  %{
  PURPOSE:
  Compute f such that f=0 for monodromy matrix

  INPUT:
 
  %}

  n = size(traj,2);
  [phi, L] = unpack_state(s);
  
  %velocity of the comoving monodromy matrix
  v = 0*phi;
  for i = 1:n
    J = jacobian( traj(:,i) );
    v(:,:,i) = J*phi(:,:,i) - phi(:,:,i)*L;
  end

  k = 0:n-1;
  k(k>n/2) = k(k>n/2) - n;
  inv_ik = 1./(1i*k);
  inv_ik(1) = 0;
  inv_ik = reshape(inv_ik, [1,1,n]);
  %No real part since phi can be complex
  integrate = @(f) ifft( inv_ik.* fft(f,[],3), [], 3 );


  f = (phi - mean(phi,3))*2*pi/T - integrate(v);
  f = reshape(f, [9*n,1]);

  init = reshape( phi(:,:,1) - eye(3), [9,1] );
  f = [f; init];
end