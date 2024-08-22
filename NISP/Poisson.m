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
%y = 1000;
y = 1e7;

%%% Omar Knio Problem
%  if( d == 1 || d == 3 || d == 4)
%       y = 1e7;
%  else
%      y = 0;    
%   end

%--- Dirichlet boundary condition ----
function y = gd(x, d, t)
y = 0;

%--- Neumann boundary condition ---
function y = gn(x, d, t)
y = 0;

%--- Right-hand side, source term ---
function y = f(x, d, t)
y =1;
% y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
%y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));
