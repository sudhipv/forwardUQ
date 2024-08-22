

%%% NON INTRUSIVE SSFEM
clearvars
%rng('default')
% Load the mesh
%square
square_refined

% Assemble matrix
A = AssembleMatrix(p, e, t, 'Poisson', [], 0);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson', [], 0);

% Solve the linear system
% U = A \ b;

%  A is being represented as a Log Normal Random Variable
n = 50000;
dim = 1;
ord_out = 10;

% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))

xi = randn(1,n);

psi(:,1) = ones(1,n);
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


A_NISP = zeros(size(A,1),size(A,2));

% Generating Realisations of Log Normal Diffusion Coefficient

         
% NISP

U = zeros(size(p,2),n_outpce);

for i = 1: n_outpce

    sum_numerator = zeros(size(p,2),1);
    
   for j = 1:n
       
       
       A_NISP = A.*(kappa(1).*psi(j,1)+kappa(2).*psi(j,2)+kappa(3).*psi(j,3)+kappa(4).*psi(j,4)+...
                    kappa(5).*psi(j,5));


                 %+...+A.*kappa(5).*psi(j,5)+A.*kappa(6).*psi(j,6) );
       b_NISP = (b.*psi(j,i));
      
      
       sum_numerator = sum_numerator + A_NISP\b_NISP;
       
       
       
   end
    
   U(:,i) = (sum_numerator/n)/factorial(i-1);
    
    
    
end

% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for ii = 2:n_outpce
    
U_std = U_std + U(:,ii).^2 .* factorial(ii-1);

end

U_std = sqrt(U_std);


% Plot solution, exact solution, and error
for j = 1:n_outpce
% % % % % figure(j); clf
% % % % % pdesurf(p,t,U(:,j))
% % % % % shading faceted
% % % % % if(j==1)
% % % % % title('Mean vector of solution expansion')
% % % % % else
% % % % %     title([num2str(j),'-PC coefficient of solution using 1 RV with 3rd order expansion',])
% % % % % end

% % % filename = sprintf('SS_NISP_%d.mat',j);
% % % pcecoeff_ssNI = U(:,j); 
% % % cd ../results
% % % save(filename,'pcecoeff_ssNI')
% % % cd ../NISP


end

% % % % % figure(n_outpce+1); clf
% % % % % pdesurf(p,t,U_std)
% % % % % shading faceted
% % % % % title('Standard deviation of solution expansion')
% % % U_std_SSNI = U_std;
% % % cd ../results
% % % save('SS_NISP_std','U_std_SSNI')
% % % cd ../NISP


% % % U_pdf(:,1) = U(9,1).*psi(:,1) + U(9,2).*psi(:,2) + U(9,3).*psi(:,3)+...
% % %              U(9,4).*psi(:,4)+ U(9,5).*psi(:,5)+ U(9,6).*psi(:,6) ;


%%% Generating PDF at centre of domain
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
U_coeff_pdf_SSNISP = zeros(n_pdf,n_outpce);

for p = 1:n_outpce
    
U_coeff_pdf_SSNISP(:,p) = (U(lo,p).*psi_pdf(:,p));   
U_pdf = U_pdf + U(lo,p).*psi_pdf(:,p);

end

figure(n_outpce+2)
[f_NI_5,xi_NI_5] = ksdensity(U_pdf);
load MCS_50000_03.mat
plot(xi_MCS,f_MCS,xi_NI_5,f_NI_5)
title('PDF of solution at centre of domain')
legend('MCS','NISP-81 nodes')


mu_U_point = mean(U_pdf)
sigma_U_point = sqrt(var(U_pdf))
skew = skewness(U_pdf)
kur = kurtosis(U_pdf)

%%% Use below code to save the values for PDF generation saved in results
%%% folder.

% % % cd ../results;
% % % save('ksd_ssfem_NI_1000_5ord.mat','xi_NI_5','f_NI_5')
% % % cd ../checks/
% % % save('PC_pdf_SSNISP.mat','U_coeff_pdf_SSNISP')
