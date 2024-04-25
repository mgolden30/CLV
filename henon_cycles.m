%{
Converge Cycles of Henon
%}

n = 7; %cycle length
trials = 100;
rng(1);

cycles = {};
i = 1;

for tr = 1:trials
  x = 1.5*(2*rand(n,1) - 1);
  [x,f] = newton(x);

  if (norm(f) < 1e-12) && is_distinct(x, cycles)
    cycles{i} = x;
    i = i+1;
  end
end

cycles

function bool = is_distinct(x, cycles)
  bool = true;
  f = @(x) [ sum(x); sum(x.^2); sum(x.^3) ];
  for i = 1:numel(cycles)
    if norm( f(x) - f(cycles{i}) ) < 1e-7
      bool = false;
      return;
    end
  end
end

function [x,f] = newton(x)
  n = numel(x);
  maxit = 1024;
  for i = 1:maxit
    f = obj(x);
    fprintf("%d: |f| = %e\n", i, norm(f) );
    if norm(f) < 1e-13
      break;
    end 

    J = zeros(n,n);
    for j = 1:n
      h = 1e-5;
      x2 =x; x2(j) = x2(j) + h;
      J(:,j) = (obj(x2) - f)/h;
    end

    [dx, ~] = lsqr( J, f );
    x = x - 0.1*dx;
  end
end



function f = obj(x)
  a = 1.4;
  b = 0.3;
  f = circshift(x,-1) - 1 + a*x.^2 - b*circshift(x,1);
end