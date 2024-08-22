

%%%% SSFEM-MCS- Poisson Equation Solver 'Kappa'/Diffusion
%%%% Coefficient modeled as a Log Normal Random Process 


%%% Generate Cijk and normsquared using Cijk.m 
%%% Generate mindex00000.dat using UQTk


clearvars
% Load the mesh
% square

% Number of Samples
n_sample = 10000;

% Generated from 'ExtractGmsh.m' 
p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';
%  square_refined

ord_in  = 2;
ord_out = 2;
dim     = 2;
str5    = '040002';
str6    = '040002';
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
end


filename_mi = strcat('../misc/cijk_ord_dim/mindex',str6,'.mat');
filename_norm = strcat('../misc/cijk_ord_dim/norm_squared',str5,'.mat');


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


U = zeros(size(p,2),n_sample);

for dd = 1:dim
 
    xi(:,dd) = randn(n_sample,1);
    
end

for nn = 1:n_sample
% Assemble matrix
A(:,:) = AssembleMatrix_ss(p, e, t, 'Poisson_ss', [], 0, s_spectral,xi(nn,:),file);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson_ss', [], 0);


% Solve the linear system
% % % % A_jk_sparse = sparse(A_jk);
% % % % [L,U] = ilu(A_jk_sparse,struct('type','ilutp','droptol',1e-5));
% % % % [U,flag,relres,iter,resvec] = bicg(A_jk,b_k,1*10^(-6),100,L,U)
%[U,flag,relres,iter,resvec] = pcg(A(:,:,1),b)


U(:,nn) = A(:,:) \ b;

end

% Plot solution, exact solution, and error
% Mean and Standard Deviation of solution
 U_mean = sum(U,2)/nn;
 
 sum_MC = zeros(size(p,2),1);
 for j = 1:nn
   
     sum_MC = sum_MC + (U(:,j) - U_mean).^2;
     
     
 end
 U_std = sum_MC/(nn-1);
 U_std = sqrt(U_std);

% Plot solution and error
figure(1); clf
pdesurf(p,t,U_mean)
shading faceted
title("Mean of Solution MC ensemble")
 
figure(2); clf
pdesurf(p,t,U_std)
shading faceted
title('Standard deviation of solution MC ensemble')

lo = 9;
[f_MCS,xi_MCS] = ksdensity(U(lo,:));
figure(4)
plot(xi_MCS,f_MCS)
title('PDF of solution at a particular point by KDE')

hold on 
hist_MCS = histogram(U(lo,:),'Normalization','pdf');
ss = sum(hist_MCS.BinWidth.*hist_MCS.Values)


mean(U(lo,:))
std(U(lo,:))

