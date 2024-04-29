function [phi, monodromy] = estimate_monodromy( zs, n, T )
  %{
  PURPOSE:
  %}

  dt  = T/n;
  phi = zeros(3,3,n+1);
  phi(:,:,1) = eye(3);
  for i = 1:n
    %index of next position
    ip  = mod(i,n)+1; %cyclic
    ip2 = i+1; %noncyclic
    
    %Use Crank-Nicolson to solve for the next timestep
    J  = jacobian( zs(:,i)  );
    Jp = jacobian( zs(:,ip) );

    phi(:,:,ip2) =         ( eye(3) + dt/2*J )*phi(:,:,i);
    phi(:,:,ip2) = linsolve( eye(3) - dt/2*Jp, phi(:,:,ip2) );
  end

  monodromy = phi(:,:,end);
  phi = phi(:,:,1:end-1);
end