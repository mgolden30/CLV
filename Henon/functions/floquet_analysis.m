function [v_e, m_e, v_c, m_c] = floquet_analysis( cycles, p )
  n = numel(cycles);

  %vectors and multipliers
  v_e = cell(n,1);
  m_e = zeros(n,1);
  v_c = cell(n,1);
  m_c = zeros(n,1);
  
  for i = 1:n
    c = cycles{i};
    c = c.'; %make size [2,p]

    %Make a mega-Jacobian
    J = zeros(2*p, 2*p);
    
    for j = 1:p
      jp = mod(j,p) + 1;
      jm = mod(j-2,p) + 1;
      
        
      b  = 2*(j -1) + (1:2); %this block
      bm = 2*(jm-1) + (1:2);
      bp = 2*(jp-1) + (1:2); 
    

      J(bp,b) = jacobian( c(:,j) );
       
    end

    [V,D] = eigs( J, 2*p );
    D = diag(D);

    %Now we have an interesting matrix. D will eigenvalues related by roots
    %of unity. Just pick the first one for each magnitude
    idx = [1; p+1];
    
    D = D(   idx);
    V = V(:, idx);  
    D
    phase = D./abs(D);

    for j = 1:p
      b  = 2*(j -1) + (1:2); %this block

      V(b,:) = V(b,:) .* phase.'.^j;
    end


    %if any( imag(V) > 0.01 )
    %  V = 1i*V;
    %end

    %Kill absolute phase!
    V = V./V(1,:);

    %Take real part
    V = real(V);

    v_e{i} = reshape( V(:,1), [2,p] );
    v_c{i} = reshape( V(:,2), [2,p] );

    D = abs(D);
    m_e = D(1);
    m_c = D(2);
  end
    

end



function v = jacobian_action( c, v )
  %c - a [2,p] cycle
  %v - a [2,p] vector field
  %Apply the Jacobian of the forward map
  
  a = 1.4; b = 0.3;
  
  %save the x component before overwriting
  vx = v(1,:); 

  v(1,:) = -2*a*c(1,:).*v(1,:) + v(2,:);
  v(2,:) = b*vx;
  
  dim = 2;
  v = circshift( v, 1, dim); %push forward
end

function v = inverse_jacobian_action( c, v )
  %c - a [2,p] cycle
  %v - a [2,p] vector field
  %Apply the Jacobian of the forward map
  
  a = 1.4; b = 0.3;
  
  %save the y component before overwriting
  vy = v(2,:); 

  v(2,:) = v(1,:) + 2*a/b/b*c(2,:);
  v(1,:) = vy/b;
  
  dim = 2;
  v = circshift( v, -1, dim); %push backwards
end