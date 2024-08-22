

%%%% SSFEM - Poisson Equation Solver 'Kappa'/Diffusion
%%%% Coefficient modeled as a Log Normal Random Process


%%% Generate Cijk and normsquared using Cijk.m 
%%% Generate mindex00000.dat using UQTk


clearvars
% Load the mesh
%  square

% Generated from 'ExtractGmsh.m'  - Don't forget to change the lines in
% Assemblematrix and AssembleVector -@ Iterating over all edges

tic

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';
% square_refined

ord_in  = 2;
ord_out = 3;
dim     = 9;

%%%% first two digits represent order of expansion and last two
%%%% digits represent random variable

% str5    = '_O3_D9';
% str6    = '_O3_D9';

str5    = '_O3_D9';
str6    = '_O3_D9';

% Inherent Gaussian Properties
mu_g = 0;
sigma_g = 0.3;

% Correlation Length
b = 1;
% Half the domain length
a = 0.5;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of input random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))

% num_spectral : determines the kind of expansion employed for input
% parameter
%%%% 1 ........ Gaussian Random Variable
%%%% 2 ........ Log Normal Random Variable
%%%% 3 ........ Gaussian Process % Not yet active
%%%% 4 ........ LogNormal Process

num_spectral = 4;

if(num_spectral >2)
    GetLambdaSymmetric_2D(b,a,dim);
    file = load('lambda_2D.mat');
else
    file = 0;
end


filename_mi = strcat('./../misc/cijk_ord_dim/mindex',str6);
filename_norm = strcat('./../misc/cijk_ord_dim/norm_squared',str5);

s_spectral = struct;
s_spectral.num_spectral = num_spectral;
s_spectral.n_inpce = n_inpce;
s_spectral.n_outpce = n_outpce;
s_spectral.dim = dim;
s_spectral.ord_in = ord_in;
s_spectral.ord_out = ord_out;
s_spectral.mu_g = mu_g;
s_spectral.sigma_g = sigma_g;
s_spectral.mindex = load(filename_mi);
s_spectral.norm = load(filename_norm);

A = zeros(size(p,2),size(p,2),n_inpce);

for ipce=1:n_inpce
% Assemble matrix
A(:,:,ipce) = AssembleMatrix_ss(p, e, t, 'Poisson_ss_Normalized', [], 0, s_spectral,ipce,file);
end

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson_ss_Normalized', [], 0);

%%% Intrusive SSFEM solve %%%%
 A_jk = zeros(size(A,1)*n_outpce);
%%% Multiplication Tensor
filename = strcat('./../misc/cijk_ord_dim/Cijk',str5);
Cijk = load(filename);
Cijk.moment_PsiNorm;
% Cijk.moment_ND;
n_A = size(A,1);

tic
 
% Assembly block struture
for k = 1: n_outpce
    for j = 1:n_outpce
        
        for i= 1: n_inpce
     
    A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A)...
        = A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A) + Cijk.moment_PsiNorm(i,j,k) * A(:,:,i);
    
        end
        
    end
    
end

toc

% Force Vector
b_k = zeros(size(b,1)*n_outpce,1);
% Expectation of psi_k <b psi_k>/<psi_k^2>

% exp_pce = [1, 0, 0, 0, 0];
exp_pce = zeros(1,n_outpce);
exp_pce(1,1) = 1;
for k = 1: n_outpce
    
    b_k(((k-1)*size(b,1)+1:(k-1)*size(b,1)+size(b,1)),1) = b(1:size(b,1),1)* exp_pce(k);
    
end

% Solve the linear system
% % % % A_jk_sparse = sparse(A_jk);
% % % % [L,U] = ilu(A_jk_sparse,struct('type','ilutp','droptol',1e-5));
% % % % [U,flag,relres,iter,resvec] = bicg(A_jk,b_k,1*10^(-6),100,L,U)
%[U,flag,relres,iter,resvec] = pcg(A(:,:,1),b)

%%% Solution Process
 U = A_jk(:,:) \ b_k;
%%% Kronecker Ajk used
% U = A_jk_Kron(:,:) \ b_k;


% Plot solution, exact solution, and error
index_coeff = size(U,1)/n_outpce;
for j = 1:n_outpce
figure(j); clf
pdesurf(p,t,U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff))
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    titlename = sprintf('%d th PC coefficient of solution using %d RV with %drd order expansion',j,dim,ord_out);
    title(titlename)
end

end


norm = s_spectral.norm;

% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for j = 2:n_outpce
U_std = U_std + U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff).^2; %.* norm.norm_squared(j)
end

U_std = sqrt(U_std);

figure(n_outpce+1); clf
pdesurf(p,t,U_std)
shading faceted
title('Standard deviation of solution Expansion')

U_CoV = U_std./U(1:index_coeff);

U_CoVmax = max(abs(U_CoV));


% % % figure(n_outpce+2); clf
% % % pdesurf(p,t,U_CoV)
% % % shading faceted
% % % title('Coefficient of Variation of solution Expansion')


%%% Plotting PDF at center of Domain

lo = 81;
n_pdf = 600000;

for dd = 1:dim
xi_pdf(:,dd) = randn(n_pdf,1);
end

psi_pdf = zeros(n_pdf,n_outpce);

                switch(dim)

                case 2
                                

                                psi_pdf(:,1) = ones(1,n_pdf);
                                psi_pdf(:,2) = xi_pdf(:,1);
                                psi_pdf(:,3) = xi_pdf(:,2);
                                psi_pdf(:,4) = xi_pdf(:,1).^2-1;
                                psi_pdf(:,5) = xi_pdf(:,1).*xi_pdf(:,2);
                                psi_pdf(:,6) = xi_pdf(:,2).^2-1;
                                psi_pdf(:,7) = xi_pdf(:,1).^3-3*xi_pdf(:,1);
                                psi_pdf(:,8) = xi_pdf(:,1).^2.*xi_pdf(:,2)-xi_pdf(:,2);
                                psi_pdf(:,9) = xi_pdf(:,1).*xi_pdf(:,2).^2-xi_pdf(:,1);
                                psi_pdf(:,10) = xi_pdf(:,2).^3-3*xi_pdf(:,2);


                case 3
                                psi_pdf(:,1) = 1.0;
                                psi_pdf(:,2) = xi_pdf(:,1);
                                psi_pdf(:,3) = xi_pdf(:,2);
                                psi_pdf(:,4) = xi_pdf(:,3);
                                psi_pdf(:,5) = (xi_pdf(:,1).^2 - 1);
                                psi_pdf(:,6) = xi_pdf(:,1).*xi_pdf(2);
                                psi_pdf(:,7) = xi_pdf(:,1).*xi_pdf(3);
                                psi_pdf(:,8) = (xi_pdf(:,2).^2 - 1);
                                psi_pdf(:,9) = (xi_pdf(:,2).*xi_pdf(:,3));
                                psi_pdf(:,10) = (xi_pdf(:,3).^2 - 1);
                                psi_pdf(:,11) = (xi_pdf(:,1).^3 - 3*xi_pdf(:,1));
                                psi_pdf(:,12) = (xi_pdf(:,1).^2.*xi_pdf(:,2) - xi_pdf(:,2));
                                psi_pdf(:,13) = (xi_pdf(:,1).^2.*xi_pdf(:,3) - xi_pdf(:,3));
                                psi_pdf(:,14) = (xi_pdf(:,2).^2.*xi_pdf(:,1) - xi_pdf(:,1));
                                psi_pdf(:,15) = xi_pdf(:,1).*xi_pdf(:,2).*xi_pdf(:,3);
                                psi_pdf(:,16) = (xi_pdf(:,3).^2.*xi_pdf(:,1) - xi_pdf(:,1));
                                psi_pdf(:,17) = (xi_pdf(2).^3 - 3*xi_pdf(:,2));
                                psi_pdf(:,18) = (xi_pdf(:,2).^2.*xi_pdf(:,3) - xi_pdf(:,3));
                                psi_pdf(:,19) = (xi_pdf(:,3).^2.*xi_pdf(:,2) - xi_pdf(:,2));
                                psi_pdf(:,20) = (xi_pdf(:,3).^3 - 3*xi_pdf(:,3));

                end

U_pdf = zeros(n_pdf,1);  

for pp = 1:n_outpce
U_pdf(:,1) = U_pdf + U(lo+(pp-1)*index_coeff,1).*psi_pdf(:,pp)/sqrt(norm.norm_squared(pp));
end

figure(n_outpce+3)
[f_IN,xi_IN] = ksdensity(U_pdf);
plot(xi_IN,f_IN)
title('PDF of solution at center of Domain')
hold on

hist_pdf = histogram(U_pdf,'Normalization','pdf');
ss = sum(hist_pdf.Values.*hist_pdf.BinWidth)

mean(U_pdf)
std(U_pdf)

