
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


% % % ydata = f(qdpts);

pc_coeff = zeros(size(p,2),ord_out+1);

for n_pc = 1:ord_out+1

    sum_f = 0;
    
    for i = 1:size(xi,1)
  
    f = ((psi(i,n_pc)*w(i))/(kappa(1)*1+ kappa(2)*xi(i,:) + kappa(3)* (xi(i,:)^2-1) +...
            kappa(4)*(xi(i,:)^3-3*xi(i,:))+kappa(5)*(xi(i,:)^4-6*xi(i,:)^2+3)));
      
        sum_f = sum_f + f;
    end
        
   
    coeff(1,n_pc) = sum_f/factorial(n_pc-1);
    pc_coeff(:,n_pc) = numerator.*coeff(1,n_pc);
   
  
end


lo = 9;
n_pdf = 300000;


xi_pdf = randn(1,n_pdf);

psi_pdf(:,1) = ones(1,n_pdf);
psi_pdf(:,2) = xi_pdf;
psi_pdf(:,3) = xi_pdf.^2-1;
psi_pdf(:,4) = xi_pdf.^3-3.*xi_pdf;
psi_pdf(:,5) = xi_pdf.^4-6*xi_pdf.^2+3;
psi_pdf(:,6) = xi_pdf.^5-10*xi_pdf.^3+15*xi_pdf;
psi_pdf(:,7) = xi_pdf.^6-15*xi_pdf.^4+45*xi_pdf.^2-15;
psi_pdf(:,8) = xi_pdf.^7-21*xi_pdf.^5+105*xi_pdf.^3-105*xi_pdf;
psi_pdf(:,9) = xi_pdf.^8-28*xi_pdf.^6+210*xi_pdf.^4-420*xi_pdf.^2+105;
psi_pdf(:,10) = xi_pdf.^9-36*xi_pdf.^7+378*xi_pdf.^5-1260*xi_pdf.^3+945*xi_pdf;
psi_pdf(:,11) = xi_pdf.^10-45*xi_pdf.^8+630*xi_pdf.^6-3150*xi_pdf.^4+4725*xi_pdf.^2-945;



U_pdf = zeros(n_pdf,1);
for pp = 1:ord_out+1
U_coeff(pp) = pc_coeff(lo,pp);
U_pdf = U_pdf + pc_coeff(lo,pp).*psi_pdf(:,pp);
end


[f_NI,xi_NI] = ksdensity(U_pdf);


% mu_U_point = mean(U_pdf)
% sigma_U_point = sqrt(var(U_pdf))
% skew = skewness(U_pdf)
% kur = kurtosis(U_pdf)

figure(200)
load ('MCS_50000_03.mat');
hold on
plot(xi_MCS,f_MCS,xi_NI,f_NI)


format long
U_coeff     
      






















% Save required data in bin folder to construct PCE surrogate, the
% filenames must be kept as is for subsequent call to pce_resp
% % % % % cd(uqtk_path);
% % % % % save ( 'ydata.dat', 'ydata', '-ascii');
% % % % % save ('xdata.dat', 'qdpts', '-ascii');
% save ('wghts.dat', 'wghts', '-ascii');

% obtain the PC coefficients for the surrogate model
% The ?-d ? option specifies the dimensionality of the problem and ?-o?
% specifies the order of PC representation. The result of executing the 
% above: a file named PCcoeff quad.dat containing the PC 
% coefficients; ydata_pc.dat conatining the PC response approximation at
% the quadrature points; mindex.dat containing the multi-indices
% % % % % % [status,out] = system(['./pce_resp -x HG -d ', num2str(dim),' -o ', num2str(ord_out)'])
% % % % % % 
% % % % % % load PCcoeff_quad.dat; %PCE coefficients
% % % % % % load mindex.dat %Multi-indices


