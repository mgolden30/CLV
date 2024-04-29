function J = jacobian(x)
  %Jacobian of the forward map
  a = 1.4; b = 0.3;
  J = [-2*a*x(1), 1; b, 0];
end