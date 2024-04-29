function J = backwards_jacobian( x )
  %{
  invert the henon map
  %}

  a= 1.4; b= 0.3;
  J = [0, 1/b;
       1, 2*a/b/b*x(2) ];
end