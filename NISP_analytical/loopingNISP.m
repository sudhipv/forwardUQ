

%%% NON INTRUSIVE SSFEM
clearvars

for num_loop = 1:1000
% square
square_refined

n = 50000;
dim = 1;
ord_out = 7;

% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim));

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


% NISP

U = zeros(size(p,2),n_outpce);

for i = 1: n_outpce

    sum_numerator = 0;
    
   for j = 1:n

       denominator = (kappa(1).*psi(j,1)+kappa(2).*psi(j,2)+kappa(3).*psi(j,3)+kappa(4).*psi(j,4)+...
                      kappa(5).*psi(j,5));
%                   +kappa(6).*psi(j,6);
%                   +kappa(7).*psi(j,7)+ kappa(8).*psi(j,8));


       NISP_Analytic = (psi(j,i))/(denominator);

          
      
      
       sum_numerator = sum_numerator + NISP_Analytic;
       
       
       
   end
   coeff = (sum_numerator/n)/factorial(i-1);
   U(:,i) = (sin(pi*p(1,:)).* sin(pi*p(2,:)))*coeff;
    
    
    
end

% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for ii = 2:n_outpce
    
U_std = U_std + U(:,ii).^2 .* factorial(ii-1);

end

U_std = sqrt(U_std);

% Plot solution, exact solution, and error
% % for j = 1:n_outpce
% % figure(j); clf
% % pdesurf(p,t,U(:,j))
% % shading faceted
% % if(j==1)
% % title('Mean vector of solution expansion')
% % else
% %     title([num2str(j),'-PC coefficient of solution using 1 RV with 3rd order expansion',])
% % end


% % % filename = sprintf('NISP_%d.mat',j);
% % % pcecoeff_NI = U(:,j); 
% % % cd ../results
% % % save(filename,'pcecoeff_NI')
% % % cd ../NISP_analytical/

% % end

U_std_NI = U_std;
% % % % cd ../results
% % % % save('NISP_std','U_std_NI')
% % % % cd ../NISP_analytical/


% % figure(n_outpce+1); clf
% % pdesurf(p,t,U_std)
% % shading faceted
% % title('Standard deviation of solution expansion')

%%%% PDF at apoint generated using KDE and histogram

lo = 9;
% % % % % % n_pdf = 300000;
% % % % % % 
% % % % % % 
% % % % % % xi_pdf = randn(1,n_pdf);
% % % % % % 
% % % % % % psi_pdf(:,1) = ones(1,n_pdf);
% % % % % % psi_pdf(:,2) = xi_pdf;
% % % % % % psi_pdf(:,3) = xi_pdf.^2-1;
% % % % % % psi_pdf(:,4) = xi_pdf.^3-3.*xi_pdf;
% % % % % % psi_pdf(:,5) = xi_pdf.^4-6*xi_pdf.^2+3;
% % % % % % psi_pdf(:,6) = xi_pdf.^5-10*xi_pdf.^3+15*xi_pdf;
% % % % % % psi_pdf(:,7) = xi_pdf.^6-15*xi_pdf.^4+45*xi_pdf.^2-15;
% % % % % % psi_pdf(:,8) = xi_pdf.^7-21*xi_pdf.^5+105*xi_pdf.^3-105*xi_pdf;
% % % % % % psi_pdf(:,9) = xi_pdf.^8-28*xi_pdf.^6+210*xi_pdf.^4-420*xi_pdf.^2+105;
% % % % % % psi_pdf(:,10) = xi_pdf.^9-36*xi_pdf.^7+378*xi_pdf.^5-1260*xi_pdf.^3+945*xi_pdf;
% % % % % % psi_pdf(:,11) = xi_pdf.^10-45*xi_pdf.^8+630*xi_pdf.^6-3150*xi_pdf.^4+4725*xi_pdf.^2-945;



% % % U_pdf = zeros(n_pdf,1);
% % % U_coeff_pdf = zeros(n_pdf,n_outpce);
for pp = 1:ord_out+1
U_coeff(pp) = U(lo,pp);
% % % U_coeff_pdf(:,pp) = (U(lo,pp).*psi_pdf(:,pp));
% % % U_pdf = U_pdf + U(lo,pp).*psi_pdf(:,pp);
end


% % % [f_NI,xi_NI] = ksdensity(U_pdf);
% % % % % plot(xi_NI,f_NI)
% % % % % title('PDF of solution at a particular point point from ensemble')


% mu_U_point = mean(U_pdf)
% sigma_U_point = sqrt(var(U_pdf))
% skew = skewness(U_pdf)
% kur = kurtosis(U_pdf)

% figure(200)
% % openfig('NISP_50000.fig');
% % hold on
% % plot(xi_NI,f_NI)


% % % figure(200+num_loop)
% % % load ('MCS_50000_03.mat');
% % % hold on
% % % plot(xi_MCS,f_MCS,xi_NI,f_NI)
%%% Use below code to save the values for PDF generation saved in results
%%% folder.

% cd ../results;
% save('analytcal2_NISP_1000_5ord.mat','xi_NI_5','f_NI_5')
format long
U_coeff;

looping_coeff_50000(num_loop,:) = U_coeff; 

end

% % % % figure(200)
% % % % h = histogram(U_pdf,1000,'Normalization','probability');
% % % % NISP_hist = h.Values;
% % % % cd ../../jupyter
% % % % save('NISP_hist.mat','NISP_hist')
% % cd ../checks/
% % save('NISP_coeff_10000.mat','looping_coeff_10000')
