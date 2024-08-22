function integral = Poisson_ss(u, v, w, du, dv, dw, dx, ds, x, d, t, eq, ord_in, dim, num_spectral,n_ipce)

% Variational formulation of Poisson's equation.
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

if eq == 1
  kappa = GenerateCoeff(ord_in, dim, num_spectral,x, d, t);  
  
  %%% Applying BC only to mean term
%     if n_ipce == 1
        integral = kappa(n_ipce) * du'*dv*dx + g(x,d,t)*u*v*ds;
%     else
%          integral = kappa(n_ipce) * du'*dv*dx;
%     end
  
else
        integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;    
  
end

% Function to generate the coefficients of Spectral Expansion for Random
% Diffusivity

function kappa = GenerateCoeff(ord_in, dim, num_spectral, x, d, t)


switch (num_spectral) 
    
    % Gaussian Random variable
    case 1 
     mu_g = 1.0;
     sigma_g = 0.3;
     kappa(1) = mu_g;
     kappa(2) = sigma_g;
     
    % Log Normal Random Variable
    case 2 
     
     mu_g = 0;
     sigma_g = 0.3;
    
     mu_l = exp(mu_g + sigma_g^2/2);

%%%% for omar knio problem
%      mu_l = 1.0;
%      sigma_g = 0.4039;

     
     kappa(1) = mu_l;
     kappa(2) = mu_l*sigma_g;
     kappa(3) = mu_l*sigma_g^2/sqrt(2);   
     kappa(4) = mu_l*sigma_g^3/sqrt(6);     
     kappa(5) = mu_l*sigma_g^4/sqrt(24);
     kappa(6) = mu_l*sigma_g^5/sqrt(120);
     kappa(7) = mu_l*sigma_g^6/sqrt(720);
     kappa(8) = mu_l*sigma_g^7/sqrt(5040);   
     kappa(9) = mu_l*sigma_g^8/sqrt(40320);
     kappa(10) = mu_l*sigma_g^9/sqrt(362880);
     kappa(11) = mu_l*sigma_g^10/sqrt(3628800);
%      kappa(3) = mu_l*sigma_g^2/2;   
%      kappa(4) = mu_l*sigma_g^3/6;     
%      kappa(5) = mu_l*sigma_g^4/24;
%      kappa(6) = mu_l*sigma_g^5/120;
%      kappa(7) = mu_l*sigma_g^6/720;
%      kappa(8) = mu_l*sigma_g^7/5040;   
        
     % Gaussian process  
     
    case 3 
      
     % Log Normal Random process
    case 4    
        
end


%--- Conductivity (penalty factor) ---
function y = g(x, d, t)
%%% d is the edge 

y = 1e15; % aLL fixed Boundary
% 1 - bottom , 2 - right , 3 - top, 4 - left
% d == 1 || d == 2 || d == 3 ||

%%% Simple plane surface solution problem
%  if( d == 2 || d == 4)
%     y = 1e7;
%  else
%      y = 0;    
% end

%%% Omar Knio Problem
%  if( d == 1 || d == 3 || d == 4)
%       y = 1e7;
%  else
%      y = 0;    
%   end





%--- Dirichlet boundary condition ----
function y = gd(x, d, t)


%%% All fixed Boundary

y =0;


% if( d == 2)
%     y = 1;
% else
%      y = 0;
% end

%%% OmAR kniO

% if( d == 1 || d == 3 || d == 4)
%     y = 0;
% else
%     y = 0;
% end






%--- Neumann boundary condition ---
function y = gn(x, d, t)
% y = 0;

%%% oMAR knio 

 y =0;


 

%--- Right-hand side, source term ---
function y = f(x, d, t)
%  y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
   y = 1;
%  y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));

