function v = rossler( x )
  a = 0.2;
  b = 0.2;
  c = 5.7;
  v = [-x(2,:)-x(3,:);
       x(1,:) + a*x(2,:);
       b + x(3,:).*(x(1,:)-c)];
end