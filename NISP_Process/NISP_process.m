

%%%% SSFEM-MCS- Poisson Equation Solver 'Kappa'/Diffusion
%%%% Coefficient modeled as a Log Normal Random Process 


%%% Generate Cijk and normsquared using Cijk.m 
%%% Generate mindex00000.dat using UQTk


clearvars
% Load the mesh
% square

% Number of quadrature points
num_qd = 5;

% Mesh from GMSH
% Generated from 'ExtractGmsh.m' 
% % % p = load('points.txt')';
% % % e =load('edges.txt')';
% % % t =load('triangles.txt')';
 square_refined

ord_in  = 2;
ord_out = 3;
dim     = 3;
str5    = '040003';
str6    = '030003';
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

filename_mi = strcat('mindex',str6,'.dat');
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

% Quadrature weights and points

% point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);


[status,out] = system(['./generate_quad -d ', num2str(dim), ' -p ', num2str(num_qd), ' -g HG -x full'])

% load the quadrature points and weights
xi = load('qdpts.dat'); % Need to scale each variable, qdpts = f(qdpts)
w = load('wghts.dat');

cd ../../../ssfem_matlab/NISP_Process

t_qd = size(xi,1);
%%% Calculating multidimensional psi values at quadrature points

psi = zeros(t_qd,n_outpce);

                switch(dim)

                case 2
                                

                                psi(:,1) = 1.0;
                                psi(:,2) = xi(:,1);
                                psi(:,3) = xi(:,2);
                                psi(:,4) = xi(:,1).^2-1;
                                psi(:,5) = xi(:,1).*xi(:,2);
                                psi(:,6) = xi(:,2).^2-1;
                                psi(:,7) = xi(:,1).^3-3*xi(:,1);
                                psi(:,8) = xi(:,1).^2.*xi(:,2)-xi(:,2);
                                psi(:,9) = xi(:,1).*xi(:,2).^2-xi(:,1);
                                psi(:,10) = xi(:,2).^3-3*xi(:,2);


                case 3
                                psi(:,1) = 1.0;
                                psi(:,2) = xi(:,1);
                                psi(:,3) = xi(:,2);
                                psi(:,4) = xi(:,3);
                                psi(:,5) = (xi(:,1).^2 - 1);
                                psi(:,6) = xi(:,1).*xi(2);
                                psi(:,7) = xi(:,1).*xi(3);
                                psi(:,8) = (xi(:,2).^2 - 1);
                                psi(:,9) = (xi(:,2).*xi(:,3));
                                psi(:,10) = (xi(:,3).^2 - 1);
                                psi(:,11) = (xi(:,1).^3 - 3*xi(:,1));
                                psi(:,12) = (xi(:,1).^2.*xi(:,2) - xi(:,2));
                                psi(:,13) = (xi(:,1).^2.*xi(:,3) - xi(:,3));
                                psi(:,14) = (xi(:,2).^2.*xi(:,1) - xi(:,1));
                                psi(:,15) = xi(:,1).*xi(:,2).*xi(:,3);
                                psi(:,16) = (xi(:,3).^2.*xi(:,1) - xi(:,1));
                                psi(:,17) = (xi(:,2).^3 - 3*xi(:,2));
                                psi(:,18) = (xi(:,2).^2.*xi(:,3) - xi(:,3));
                                psi(:,19) = (xi(:,3).^2.*xi(:,2) - xi(:,2));
                                psi(:,20) = (xi(:,3).^3 - 3*xi(:,3));

                end


norm = s_spectral.norm;


U = zeros(size(p,2),n_outpce);

for n_pc = 1:n_outpce
    
sum_numerator = zeros(size(p,2),1);

for np = 1:t_qd
    
% Assemble matrix
A(:,:) = AssembleMatrix_ss(p, e, t, 'Poisson_ss', [], 0, s_spectral,xi(np,:),file);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson_ss', [], 0);

 b = (b.*psi(np,n_pc)*w(np));
    
 sum_numerator = sum_numerator + A\b;


end

U(:,n_pc) = (sum_numerator)/sqrt(norm.norm_squared(n_pc));

end


for j = 1:n_outpce
figure(j); clf
pdesurf(p,t,U(:,j))
shading interp
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    titlename = sprintf('%d th PC coefficient of solution using %d RV with %drd order expansion',j,dim,ord_out);
    title(titlename)
end

end


% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for ii = 2:n_outpce
    
U_std = U_std + U(:,ii).^2;

end

U_std = sqrt(U_std);

figure(n_outpce+1); clf
pdesurf(p,t,U_std)
shading faceted
title('Standard deviation of solution Expansion')

U_CoV = U_std./U(:,1);

U_CoVmax = max(abs(U_CoV));




%%% Plotting PDF at center of Domain

lo = 9;
n_pdf = 300000;

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
                                psi_pdf(:,17) = (xi_pdf(:,2).^3 - 3*xi_pdf(:,2));
                                psi_pdf(:,18) = (xi_pdf(:,2).^2.*xi_pdf(:,3) - xi_pdf(:,3));
                                psi_pdf(:,19) = (xi_pdf(:,3).^2.*xi_pdf(:,2) - xi_pdf(:,2));
                                psi_pdf(:,20) = (xi_pdf(:,3).^3 - 3*xi_pdf(:,3));

                end

U_pdf = zeros(n_pdf,1);  

for pp = 1:n_outpce
U_pdf(:,1) = U_pdf + U(lo,pp).*psi_pdf(:,pp)/sqrt(norm.sorm_squared(pp));
end

figure(n_outpce+3)
[f_IN,xi_IN] = ksdensity(U_pdf);
plot(xi_IN,f_IN)
title('PDF of solution at center of Domain')
% % hold on
% % 
% % hist_pdf = histogram(U_pdf,'Normalization','pdf');
% % ss = sum(hist_pdf.Values.*hist_pdf.BinWidth)

mean(U_pdf)
std(U_pdf)




















