function x = generate_backwards_traj( x0, n )
  %{
  invert the henon map
  %}

  x = zeros(2,n+1);
  x(:,1) = x0;

  for i = 1:n
    x(:,i+1) = inverse_map(x(:,i));
  end
end

function x = inverse_map(x)
  %Inverse of henon
  a= 1.4; b= 0.3;

  temp = x(1,:); %save x
  x(1,:) = x(2,:)/b;
  x(2,:) = temp - 1 + a/b/b*x(2,:).^2;
end