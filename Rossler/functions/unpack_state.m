function [phi, L] = unpack_state( s )
  n = (numel(s)/9 - 1);

  phi = reshape( s(1:3*3*n),       [3,3,n] );
  L   = reshape( s(3*3*n+1:end)  , [3,3]   );
end