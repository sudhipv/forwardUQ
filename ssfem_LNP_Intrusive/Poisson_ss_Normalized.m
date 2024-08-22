function integral = Poisson_ss_Normalized(u, v, w, du, dv, dw, dx, ds, x, d, t, eq, s_spectral,ipce,file)

% Variational formulation of Poisson's equation.
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

if eq == 1
  kappa = GenerateCoeff(s_spectral,x, d, t,file,ipce); 
  
  %%% All PC terms having BC.
%   if ipce == 1
    integral = kappa * du'*dv*dx + g(x,d,t)*u*v*ds;
%   else
%       integral = kappa * du'*dv*dx;
%   end
  


else
  integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end

% Function to generate the coefficients of Spectral Expansion for Random
% Diffusivity

function kappa = GenerateCoeff(s_spectral, x, d, t,file,ipce)


switch (s_spectral.num_spectral) 
    
    % Gaussian Random variable
    case 1 
     
     kappa = s_spectral.mu_g;
     kappa = s_spectral.sigma_g;
     
    % Log Normal Random Variable
    case 2 
     
     mu_g = s_spectral.mu_g;
     sigma_g = s_spectral.sigma_g;
     mu_l = exp(mu_g + sigma_g^2/2);
     cd(1) = mu_l;
     cd(2) = mu_l*sigma_g;
     cd(3) = mu_l*sigma_g^2/sqrt(2);   
     cd(4) = mu_l*sigma_g^3/sqrt(6);  
     
%       kappa = cd(ipce);
     
     kappa = 1;
     
%      kappa(3) = mu_l*sigma_g^2/2;   
%      kappa(4) = mu_l*sigma_g^3/6;     
        
        
     % Gaussian process  
     
    case 3 
      
     % Log Normal Random process
    case 4    
        
        
        [kappa] = GetKLEterms_Normalized(x,d,t,file,s_spectral,ipce);
        
% % %     cfolder = pwd;
% % %         % % % point to UQTk directopry
% % %     uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
% % %     cd(uqtk_path);
% % %         
% % %     [status,out] = system(['./gen_mi -p', num2str(s_spectral.ord_out), ' -q', num2str(s_spectral.dim) , '-x TO'])

        
end


%--- Conductivity (penalty factor) ---
function y = g(x, d, t)
y = 1e7;
%y = 0.001;

%--- Dirichlet boundary condition ----
function y = gd(x, d, t)
y = 0;

%--- Neumann boundary condition ---
function y = gn(x, d, t)
y = 0;

%--- Right-hand side, source term ---
function y = f(x, d, t)
 y = 1;
 %y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
%  y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));