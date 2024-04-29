function x = generate_traj( x0, n )
  %{
  invert the henon map
  %}

  x = zeros(2,n+1);
  x(:,1) = x0;

  for i = 1:n
    x(:,i+1) = map(x(:,i));
  end
end

function x = map(x)
  %Inverse of henon
  a= 1.4; b= 0.3;

  temp = x(1,:); %save x
  x(1,:) = 1 - a*x(1,:).^2 + x(2,:);
  x(2,:) = b*temp;
end