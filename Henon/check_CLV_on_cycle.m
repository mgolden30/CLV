%{
Apply the CLV routine to check if we get the same result as from Floquet
analysis
%}

clear;

p = 15; %cycle length
cycles = lift_cycles(p);

%Traditional Floquet analysis
maxit = 1024;
[v_e, m_e, v_c, m_c] = floquet_analysis( cycles, p );

%Let's just look at the first cycle of length p
c = cycles{1};
c = c.';

v_c   = v_c{1}; %Take the Floquet analysis answer
CLV_c = v_c*0;

v_e   = v_e{1}; %Take the Floquet analysis answer
CLV_e = v_e*0;


for m = 1:p
  x0  = c(:,m);
  n     = 8;  %time forward
  maxit = 4; %power iterations 

  %Compute the orthonormal vectors
  [vv_c, vv_e] = compute_CLVs( x0, n, maxit );

  CLV_c(:,m) = vv_c;
  CLV_e(:,m) = vv_e;
end


v_c = v_c ./ vecnorm( v_c );
v_e = v_e ./ vecnorm( v_e );

%Check if this is proportional to the Floquet vector
ratio_c = CLV_c ./ v_c
ratio_e = CLV_e ./ v_e