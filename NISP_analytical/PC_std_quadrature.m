

%%%% 


clearvars
%%%% NISP using Quadrature
square_refined

ord_out = 7;
dim = 1;
num_qd = 20;

mu_g = 1;
sigma_g = 0.3;
mu_l = exp(mu_g + sigma_g^2/2);
kappa(1) = mu_l;
kappa(2) = mu_l*sigma_g;
kappa(3) = mu_l*sigma_g^2/2;   
kappa(4) = mu_l*sigma_g^3/6;
kappa(5) = mu_l*sigma_g^4/24;
kappa(6) = mu_l*sigma_g^5/120;
kappa(7) = mu_l*sigma_g^6/720;
kappa(8) = mu_l*sigma_g^7/5040;
kappa(9) = mu_l*sigma_g^8/40320;
kappa(10) = mu_l*sigma_g^9/362880;
kappa(11) = mu_l*sigma_g^10/3628800;



% Quadrature weights and points

% point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);


[status,out] = system(['./generate_quad -d ', num2str(dim), ' -p ', num2str(num_qd), ' -g HG -x full'])

% load the quadrature points and weights
xi = load('qdpts.dat'); % Need to scale each variable, qdpts = f(qdpts)
w = load('wghts.dat');

cd ../../../ssfem_matlab/NISP_analytical

numerator = (sin(pi*p(1,:)).* sin(pi*p(2,:)));
numerator_std = (sin(pi*p(1,:)).* sin(pi*p(2,:))).^2;


psi(:,1) = ones(1,size(xi,1));
psi(:,2) = xi;
psi(:,3) = xi.^2-1;
psi(:,4) = xi.^3-3.*xi;
psi(:,5) = xi.^4-6*xi.^2+3;
psi(:,6) = xi.^5-10*xi.^3+15*xi;
psi(:,7) = xi.^6-15*xi.^4+45*xi.^2-15;
psi(:,8) = xi.^7-21*xi.^5+105*xi.^3-105*xi;
psi(:,9) = xi.^8-28*xi.^6+210*xi.^4-420*xi.^2+105;
psi(:,10) = xi.^9-36*xi.^7+378*xi.^5-1260*xi.^3+945*xi;
psi(:,11) = xi.^10-45*xi.^8+630*xi.^6-3150*xi.^4+4725*xi.^2-945;


% (sin(pi*x)*sin(pi*y))/c(theta) using intrusive expansion of lognormal random variable

% %  intrsv_coeff = [0.383572209523979,-0.115066464903918,0.0172394042019312,...
% %     -0.00168384424748072,8.32383488314021e-05,2.01287144828188e-05,...
% %     -8.44655221201909e-06,1.38199882269548e-06];


pc_coeff_std = zeros(size(p,2),ord_out+1);

load PC_coeff_quad.mat

for n_pc = 1:ord_out+1

    sum_f_std = 0;
    sum_f = 0;
    
    for i = 1:size(xi,1)
  
      var_coeff = (psi(i,n_pc)/(kappa(1)*1+ kappa(2)*xi(i,:) + kappa(3)* (xi(i,:)^2-1) +...
            kappa(4)*(xi(i,:)^3-3*xi(i,:))+kappa(5)*(xi(i,:)^4-6*xi(i,:)^2+3)))^2; 
        
      f_std = (w(i)*var_coeff);
      sum_f_std = sum_f_std + f_std; 
      
     
   
    end
    
    coeff_std(1,n_pc) = sum_f_std/factorial(n_pc-1)^2;
    
    
    pc_coeff_std(:,n_pc) = numerator_std.*coeff_std(1,n_pc);
     
     
     pc_coeff_std(:,n_pc) = pc_coeff_std(:,n_pc) - pc_coeff(:,n_pc).^2;
     
     pc_coeff_std(:,n_pc) = sqrt(abs(pc_coeff_std(:,n_pc)));
     
  
end

 
 save('PC_coeff_std.mat','pc_coeff_std')
 

