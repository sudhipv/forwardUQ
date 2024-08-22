function integral = Poisson(u, v, w, du, dv, dw, dx, ds, x, d, t, eq)

% Variational formulation of Poisson's equation.
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

if eq == 1
  integral = du'*dv*dx + g(x,d,t)*u*v*ds;
else
  integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end

%--- Conductivity (penalty factor) ---
function y = g(x, d, t)
%%%% ds is zero inside boundary and takes value only at the edges
%%% d is the edge 

% 1 - bottom , 2 - right , 3 - top, 4 - left
% d == 1 || d == 2 || d == 3 ||
 if( d == 2 || d == 4)
     y = 1e7;
 else
     y = 0;    
  end

%--- Dirichlet boundary condition ----
function y = gd(x, d, t)


%%% Plane surface problem
if( d == 2)
    y = 1;
else
    y = 0;
end

%--- Neumann boundary condition ---
function y = gn(x, d, t)
y = 0;

%--- Right-hand side, source term ---
function y = f(x, d, t)
% y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
y = 1;
