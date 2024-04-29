function s = pack_state( phi, L )
  phi = reshape( phi, [numel(phi),1] );
  L   = reshape( L  , [numel(L),  1] );
  
  s = [phi; L];
end