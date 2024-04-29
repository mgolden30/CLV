function f = variational_objective( zs )
  %{
  PURPOSE:
  Compute f such that f=0 for periodic orbits

  INPUT:
  zs - 3n+1 vector (trajectory and period)
  %}

  n = (numel(zs)-1)/3;
  
  T  = zs(end);
  zs = reshape( zs(1:end-1), [3, n] );

  v = rossler(zs);

  k = 0:n-1;
  k(k>n/2) = k(k>n/2) - n;
  
  inv_ik = 1./(1i*k);
  inv_ik(1) = 0;

  integrate = @(f) real(ifft( inv_ik.* fft(f,[],2), [], 2 ));


  f = (zs - mean(zs,2))*2*pi/T - integrate(v);
  f = reshape(f, [3*n,1]);

  %Add a phase condition
  phase = zs(1,1);
  f = [f; phase];
end