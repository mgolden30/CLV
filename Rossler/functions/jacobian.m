function J = jacobian( x )
  a = 0.2;
  b = 0.2;
  c = 5.7;
  J = [0, -1, -1;
       1, a, 0;
       x(3), 0, x(1)-c];
end