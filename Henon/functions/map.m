function x = map(x)
  a= 1.4; b= 0.3;
  temp = x(1,:);
  x(1,:) = 1 - a*x(1,:).*x(1,:) + x(2,:);
  x(2,:) = b*temp;
end