function xs = rk4( x, dt, n )
  xs = zeros(3,n+1);
  xs(:,1) = x;
  for i=1:n

    k1 = dt*rossler(x);
    k2 = dt*rossler(x + k1/2);
    k3 = dt*rossler(x + k2/2);
    k4 = dt*rossler(x + k3);

    x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    xs(:,i+1) = x;
  end
  
end