
%%% Intrusive and NISP not matching after 5th coefficient even with sigma =
%%% 0.1 and 4th input and 10 th order output. Check norm.normproduct again.
%%% Howver, both closely match with OMar knio up to 4th pc coefficient.



%%%% SSFEM - Poisson Equation Solver with Random Variable 'Kappa'/Diffusion
%%%% Coefficient

%%%%NOTES

% 1. Changing Cijk to be not normalized  ?? Normalized used.


clearvars
% Load the mesh
%  square
% Generated from 'ExtractGmsh.m'  - Don't forget to change the lines in
% Assemblematrix and AssembleVector -@ Iterating over all edges

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';

% square_refined

ord_in  = 2;
ord_out = 3;
dim     = 1;
% str5    = '0100001';

str5    = '_O7_D1';

%%% for omar knio problem
% str5    = '_O4_D1';


% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))

U = zeros(size(p,2),1);

% num_spectral : determines the kind of expansion employed for input
% parameter

%%%% 2 ........ Log Normal Random Variable
%

num_spectral = 2;

for ipce=1:n_inpce
% Assemble matrix
A(:,:,ipce) = AssembleMatrix_ss(p, e, t, 'Poisson_ss', [], 0, ord_in, dim, num_spectral,ipce);
end

b = AssembleVector(p, e, t, 'Poisson_ss', [], 0);

%%% Intrusive SSFEM solve %%%%
A_jk = zeros(size(A,1)*n_outpce);
%%% Multiplication Tensor
filename = strcat('./../misc/cijk_ord_dim/Cijk',str5);

% filename = strcat('./../misc/cijk_HPC/cijk',str5);
Cijk = load(filename);
% Cijk.moment_PsiNorm;
n_A = size(A,1);

% Assembly block struture
for k = 1: n_outpce
    for j = 1:n_outpce
        
        for i= 1: n_inpce
     
    A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A)...
        = A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A) + Cijk.moment_PsiNorm(i,j,k) * A(:,:,i);
    
        end
        
    end
    
end

% Force Vector
b_k = zeros(size(b,1)*n_outpce,1);
% Expectation of psi_k <b psi_k>/<psi_k^2>
% Norm Squared

cd /Users/sudhipv/documents/ssfem_matlab/Intrusive
%%% Multiplication Tensor
%%% This is only for one dimensional expansion, if the dimension is more
%%% than 1 re calculate this norm '
norm = load('./norm_squared_O11_D1.mat');
norm.norm_squared;

% exp_pce = [1, 0, 0, 0, 0, 0,0,0];
% for k = 1: n_outpce
%     
%     b_k(((k-1)*size(b,1)+1:(k-1)*size(b,1)+size(b,1)),1) = b(1:size(b,1),1)* exp_pce(k);
%     
% end


b_k(1:size(b,1),1) = b(1:size(b,1),1);

% Solve the linear system
% % % % A_jk_sparse = sparse(A_jk);
% % % % [L,U] = ilu(A_jk_sparse,struct('type','ilutp','droptol',1e-5));
% % % % [U,flag,relres,iter,resvec] = bicg(A_jk,b_k,1*10^(-6),100,L,U)
%[U,flag,relres,iter,resvec] = pcg(A(:,:,1),b)
U = A_jk(:,:) \ b_k;



% Plot solution, exact solution, and error
index_coeff = size(U,1)/n_outpce;
for j = 1:n_outpce
figure(j+100); clf
pdesurf(p,t,U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff))
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    title([num2str(j),'-PC coefficient of solution using 1 RV with 7th order output expansion',])
end

%%% Saving PC coefficients to results folder 
% % % % % filename = sprintf('ssfem_intrsv_%d.mat',j);
% % % % % pce_coeff = U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff); 
% % % % % cd ../results
% % % % % save(filename,'pce_coeff')
% % % % % cd ../Intrusive/

end

% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for j = 2:n_outpce
U_std = U_std + U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff).^2;
end

U_std = sqrt(U_std);
figure(n_outpce+200); clf
pdesurf(p,t,U_std)
shading faceted
title('Standard deviation of solution Expansion')
ssfem_intrsv_std = U_std;
% % % cd ../results
% % % save('ssfem_intrsv_std','ssfem_intrsv_std')
% % % cd ../Intrusive/

%%

lo = 81;
n_pdf = 300000;
xi_pdf = randn(n_pdf,1);

psi_pdf(:,1) = ones(1,n_pdf);
psi_pdf(:,2) = xi_pdf;
psi_pdf(:,3) = (xi_pdf.^2-1);
psi_pdf(:,4) = (xi_pdf.^3-3.*xi_pdf);
psi_pdf(:,5) = (xi_pdf.^4-6*xi_pdf.^2+3);
psi_pdf(:,6) = (xi_pdf.^5-10*xi_pdf.^3+15*xi_pdf);
psi_pdf(:,7) = (xi_pdf.^6-15*xi_pdf.^4+45*xi_pdf.^2-15);
psi_pdf(:,8) = (xi_pdf.^7-21*xi_pdf.^5+105*xi_pdf.^3-105*xi_pdf);
psi_pdf(:,9) = xi_pdf.^8-28*xi_pdf.^6+210*xi_pdf.^4-420*xi_pdf.^2+105;
psi_pdf(:,10) = xi_pdf.^9-36*xi_pdf.^7+378*xi_pdf.^5-1260*xi_pdf.^3+945*xi_pdf;
psi_pdf(:,11) = xi_pdf.^10-45*xi_pdf.^8+630*xi_pdf.^6-3150*xi_pdf.^4+4725*xi_pdf.^2-945;




% % % %  U_pdf_og(:,1) = U(lo,:).*ones(n_pdf,1) + U(lo+index_coeff,:).*xi_pdf + U(lo+2*index_coeff,:).*(xi_pdf.^2-1)+...
% % % %              U(lo+3*index_coeff,:).*(xi_pdf.^3-3*xi_pdf)+U(lo+4*index_coeff,:).*(xi_pdf.^4-6*xi_pdf.^2+3)+...
% % % %              U(lo+5*index_coeff,:).*(xi_pdf.^5-10*xi_pdf.^3+15*xi_pdf)+...
% % % %              U(lo+6*index_coeff,:).*(xi_pdf.^6-15*xi_pdf.^4+45*xi_pdf.^2-15)+...
% % % %              U(lo+7*index_coeff,:).*(xi_pdf.^7-21*xi_pdf.^5+105*xi_pdf.^3-105*xi_pdf);
% % % %          
 
% norm = load('norm_squared070001.mat', 'norm_squared');

U_pdf = zeros(n_pdf,1);  
U_coeff_pdf_SSIN = zeros(n_pdf,n_outpce);
for pp = 1:n_outpce
U_coeff(pp) = U(lo+(pp-1)*index_coeff,1);    
U_coeff_pdf_SSIN(:,pp) = (U(lo+(pp-1)*index_coeff,1).*psi_pdf(:,pp));   
U_pdf(:,1) = U_pdf + U(lo+(pp-1)*index_coeff,1).*psi_pdf(:,pp)/sqrt(norm.norm_squared(pp));
end


figure(n_outpce+300)
[f_IN,xi_IN] = ksdensity(U_pdf);
% load 'MCS_50000_03.mat'
plot(xi_IN,f_IN)
title('PDF of solution at a point of domain by KDE')
legend('SSFEM -Intrusive')

%%% Moments of PDF
mu_U_point = mean(U_pdf)
sigma_U_point = sqrt(var(U_pdf))
skew = skewness(U_pdf)
kur = kurtosis(U_pdf)


%%% Use below code to save the values for PDF generation saved in results
%%% folder.

% % % cd ../results;
% % % save('ksd_ssfem_IN_3ord.mat','xi_IN','f_IN')

% % % cd ../checks/
% % % save('PC_pdf_SSIN.mat','U_coeff_pdf_SSIN')

% % format long
% % U_coeff
