function cycles_2d = lift_cycles( p )
  %{
  PURPOSE:
  Henon cycles are converged and stored in only the x coordinate. The y
  coordinate can be recovered.
  %}

  %Load the cell array "cycles" of desired length
  load("prime_cycles/" + p + ".mat");
  
  b = 0.3;
  n = numel( cycles );
  cycles_2d = cell( n, 1 );
  
  for i = 1:n
    c = cycles{i};
 
    %restore the y component
    c2 = [ c, b*circshift(c,1,1) ];

    cycles_2d{i} = c2; 
  end
end