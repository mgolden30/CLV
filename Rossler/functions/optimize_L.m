function s = optimize_L( s, traj, T )
  %{
  PURPOSE:
  Compute the exact best L
  %}

  n = size(traj,2);
  [phi, ~] = unpack_state(s);
  
  %velocity of the comoving monodromy matrix
  v = 0*phi;
  for i = 1:n
    J = jacobian( traj(:,i) );
    v(:,:,i) = J*phi(:,:,i);
  end

  k = 0:n-1;
  k(k>n/2) = k(k>n/2) - n;
  inv_ik = 1./(1i*k);
  inv_ik(1) = 0;
  inv_ik = reshape(inv_ik, [1,1,n]);
  %No real part since phi can be complex
  integrate = @(f) ifft( inv_ik.* fft(f,[],3), [], 3 );


  f = (phi - mean(phi,3))*2*pi/T - integrate(v);
  p = integrate(phi);

  %So L should satisfy f = p*L 
  f = reshape(f, [3*n,3]);
  p = reshape(p, [3*n,3]);
  
  %Find optimial L
  L = p\f;

  s = pack_state(phi, L);
end