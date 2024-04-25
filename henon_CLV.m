%{
CLV attempt
%}


%% Look at attractor
N = 2^13;
xs = zeros(2,N);
xs(:,1) = [0.5;0.2;];
for i = 2:N
  xs(:,i) = map( xs(:,i-1) );
end

ms = 10; %marker size
scatter( xs(1,:), xs(2,:), ms, 'filled', 'MarkerFaceColor', 'black' );

xlim([-1,1]*1.5);
ylim([-1,1]*0.5);
pbaspect([3,1,1])
axis square

%% Find some periodic orbits, since we understand these pretty well
seed = 1;
rng(seed);

cycle_length = 31; %something prime is good.

%Pick any random subset of our chaotic trajectory
start = randi(N - cycle_length);
z = xs(:,start + (1:cycle_length) - 1);
z = reshape(z, [numel(z),1]);

maxit = 1024;
h = 1e-6;
hook = 0.02;
for i = 1:maxit
  f = po_obj(z);
  J = zeros(numel(f),numel(z));
  for j = 1:numel(z)
    z2 = z; z2(j) = z2(j) + h;
    J(:,j) = (po_obj(z2) - f)/h;
  end
  [dz, ~] = lsqr( J, f );
  if( norm(f) < 1e-6)
    hook = 1;
  end
  z = z - hook*dz;
  fprintf("%d: |f| = %e\n", i, norm(f) );
end


hold on
z2 = reshape(z, [2, cycle_length]);
scatter( z2(1,:), z2(2,:), 2*ms, "MarkerFaceColor", "red" );
hold off


%% Compute the Floquet vectors and multipliers for this long orbit

%Do power iteration
v_expand = 2*rand(2,cycle_length)-1;

tic
maxit = 1024;
for i = 1:maxit
  %apply the Jacobian at each point
  for j = 1:cycle_length
    v_expand(:,j) = jacobian(z2(:,j)) * v_expand(:,j);
  end
  %cycle to next point
  v_expand = circshift(v_expand, -1, 2);
  v_expand = v_expand / max(max(abs(v_expand)));
end
toc

hold on
%make them all unit vectors for plotting!
v_expand = v_expand ./ vecnorm(v_expand);
scale= 0.05;
nice_quiver( z2(1,:), z2(2,:),  v_expand(1,:),  v_expand(2,:), scale, lw, "red" );
nice_quiver( z2(1,:), z2(2,:), -v_expand(1,:), -v_expand(2,:), scale, lw, "red" );
hold off

%Do the same thing for contracting direction
v_contract = 2*rand(2,cycle_length)-1;
tic
for i = 1:maxit
  %apply the Jacobian at each point
  for j = 1:cycle_length
    v_contract(:,j) = inv(jacobian(z2(:,j))) * v_contract(:,j);
  end
  %cycle to next point
  v_contract = circshift(v_contract, 1, 2); %cycle the other way since we are going backwards time
  v_contract = v_contract / max(max(abs(v_contract)));
end
toc

hold on
%make them all unit vectors for plotting!
v_contract = v_contract ./ vecnorm(v_contract);
lw = 2;
%scale = 0.1;
nice_quiver( z2(1,:), z2(2,:),  v_contract(1,:),  v_contract(2,:), scale, lw, "blue" );
nice_quiver( z2(1,:), z2(2,:), -v_contract(1,:), -v_contract(2,:), scale, lw, "blue" );
hold off



%% Test 
clv_contract = 0*z2;
for m = 1:cycle_length
  x0 = z2(:,m);
  n  = 128; %function applications
  V  = graham_schmidt_vectors( x0, n );

  clv_contract(:,m) = V(:,2);
end

%WHY DO I NEED THIS?!?!?!?!
%clv_contract = circshift(clv_contract, 1, 2);

%Check if this is proportional to the Floquet vector
ratio = clv_contract ./ v_contract
%SHould be all 1s and -1s!

%% Since test worked, I think I am ready to plot as a vector field

gpx = 128*4;
gpy = 128*4;
xg = linspace(-1.5,1.5,gpx);
yg = linspace(-0.5,0.5,gpy);
[xx,yy] = meshgrid(xg,yg);

clv_contract = zeros(gpy,gpx,2);

%function evaluations in the future
n = 64;
for i = 1:gpx
  for j = 1:gpy
    x0 = [xg(i); yg(j)];
    tic
    V = graham_schmidt_vectors( x0, n );
    toc
    clv_contract(j,i,:) = V(:,2);
  end
end

hold on
gray = [0.5 0.5 0.5];
lw   = 1.5;
scale = 0.25;
nice_quiver(xx, yy, clv_contract(:,:,1), clv_contract(:,:,2), scale, lw, gray);
nice_quiver(xx, yy,-clv_contract(:,:,1),-clv_contract(:,:,2), scale, lw, gray);
hold off

xlabel("x");
ylabel("y");

xticks([-1.5, 0, 1.5]);
yticks((-1:1)*0.5);

legend({"", "", "expanding floquet", "", "contracting floquet", "", "contracting CLV"})

%%
%{
openfig("my_fig.fig");

x = linspace(-1,1, 128);

hold on
  x2 = x + 0.06;
  y  = (abs(x2)).^2 + (abs(x2)).^4 - 0.12;
  %plot( x,y, "linewidth", 2 )

  %w = [-1.098; 0.00393701];
  w = [x;y];
  for i = 1:2
  scatter( w(1,:), w(2,:), "filled" );
  w = map(w);
  end
  hold off
%}

%%

n1 = clv_contract(:,:,1);
n2 = clv_contract(:,:,2);

Q11 = 1/2*(n1.^2 - n2.^2);
Q12 = n1.*n2;

sigma = 1;
Q11 = imgaussfilt(Q11, sigma);
Q12 = imgaussfilt(Q12, sigma);


drawnow;
S = Q11.^2 + Q12.^2;
S( ~isfinite(S) ) = 10000;
imagesc(xg,yg, S);
clim([0, 0.2]);
c = colorbar();
set(c, "xtick", [0:0.1:0.2]);
xticks([-1:1]*1.5);
yticks([-1:1]*0.5);

set(gca, "ydir", "normal");
colormap parula;

hold on
scatter( xs(1,:), xs(2,:), 3,'filled', 'MarkerFaceColor', 'black');
hold off
axis square

xlabel("x");
ylabel("y");
title("nematic scalar order parameter");

function V = graham_schmidt_vectors( x0, n )
  %{
  PURPOSE:
  Find the singular vectors of J with power iteration, where J is the
  jacobian of f^n.
  %} 
  
  %First generate a trajectory forward
  xs = zeros(2,n+1);
  xs(:,1) = x0;
  for i = 1:n
    xs(:,i+1) = map( xs(:,i) );
  end

  %Start with an initially orthogonal
  V = eye(2);
  maxit = 8; %I doubt I need that many iterations.
  for i = 1:maxit
    for j = 1:n
      jac = jacobian(xs(:,j));
      V   = jac * V;
      %orthonormalize
      V(:,1) = V(:,1) / vecnorm(V(:,1));
      V(:,2) = V(:,2) - V(:,1)*dot(V(:,1),V(:,2));
      V(:,2) = V(:,2) / vecnorm(V(:,2));
    end

    for j = n:-1:2
      V = jacobian(xs(:,j))' * V; %transpose backwards in time
      %orthonormalize
      V(:,1) = V(:,1) / vecnorm(V(:,1));
      V(:,2) = V(:,2) - V(:,1)*dot(V(:,1),V(:,2));
      V(:,2) = V(:,2) / vecnorm(V(:,2));
    end
  end
end


function J = jacobian(x)
  %Jacobian of the forward map
  a = 1.4; b = 0.3;
  J = [-2*a*x(1), 1; b, 0];
end

function f = po_obj(z)
  z  = reshape(z, [2, numel(z)/2]);
  zp = map(z);
  f  = circshift(z,1,2) - zp;
  f  = reshape(f, [numel(f),1]);
end

function x = map(x)
  a= 1.4; b= 0.3;
  temp = x(1,:);
  x(1,:) = 1 - a*x(1,:).*x(1,:) + x(2,:);
  x(2,:) = b*temp;
end

function nice_quiver(x,y,vx,vy,scale,lw,color)
  flatten = @(x) reshape(x, [numel(x),1]);

  x  = flatten(x);
  y  = flatten(y);
  vx = flatten(vx);
  vy = flatten(vy);

  normv = sqrt(vx.^2 + vy.^2);
  vx = vx./normv;
  vy = vy./normv;

  %{
  xl = xlim;
  yl = ylim;
  r = (yl(2)-yl(1))/(xl(2)-xl(1));
  vx = vx/r; %rescale by axis ratio to make vectors looks uniform in size
  %}

  %scale = 0.3;
  quiver(x, y, vx, vy, scale, "linewidth", lw, "color", color);
end