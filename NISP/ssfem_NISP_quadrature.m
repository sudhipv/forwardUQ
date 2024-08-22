
%%% Values don't match intrusive for higher order PCs
%%% NON INTRUSIVE SSFEM using quadrature
clearvars
% square


% Generated from 'ExtractGmsh.m'  - Don't forget to change the lines in
% Assemblematrix and AssembleVector -@ Iterating over all edges

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';

% square_refined

% Assemble matrix
A = AssembleMatrix(p, e, t, 'Poisson', [], 0);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson', [], 0);

% Solve the linear system
% U = A \ b;

%  A is being represented as a Log Normal Random Variable
dim = 1;
ord_out = 4;
n_outpce = ord_out +1;
num_qd = 100;


mu_g = 0;
% sigma_g = 0.1;
% mu_l = exp(mu_g + sigma_g^2/2);

%%%% for omar knio problem
     mu_l = 1.0;
     sigma_g = 0.4039;

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


% Quadrature weights and points

% point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);


[status,out] = system(['./generate_quad -d ', num2str(dim), ' -p ', num2str(num_qd), ' -g HG -x full'])

% load the quadrature points and weights
xi = load('qdpts.dat'); % Need to scale each variable, qdpts = f(qdpts)
w = load('wghts.dat');

cd ../../../ssfem_matlab/NISP


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


U = zeros(size(p,2),ord_out+1);

norm = load('norm_squared_O11_D1.mat', 'norm_squared');

for n_pc = 1:ord_out+1

    sum_numerator = zeros(size(p,2),1);
    
    for i = 1:size(xi,1)
  
    A_NISP = A.*(kappa(1)*1+ kappa(2)*xi(i,:) + kappa(3)* (xi(i,:)^2-1) +...
            kappa(4)*(xi(i,:)^3-3*xi(i,:))+kappa(5)*(xi(i,:)^4-6*xi(i,:)^2+3));
      
    b_NISP = (b.*psi(i,n_pc)*w(i));
    
    sum_numerator = sum_numerator + A_NISP\b_NISP;
    
    end
        
    U(:,n_pc) = (sum_numerator)/norm.norm_squared(n_pc);
  
end

% Plot solution, exact solution, and error
for j = 1:n_outpce
figure(j); clf
pdesurf(p,t,U(:,j))
shading faceted
if(j==1)
title('Mean vector -nisp')
else
    titlename = sprintf('%d - PC coefficient of solution using %d RV with %d order expansion-nisp',j,dim, ord_out);
    title(titlename)
end
end



% Calculation of standard deviation
U_std = zeros(size(p,2),1);
for ii = 2:n_outpce
    
U_std = U_std + U(:,ii).^2 .* norm.norm_squared(ii);

end

U_std = sqrt(U_std);


figure(n_outpce+1); clf
pdesurf(p,t,U_std)
shading faceted
title('Standard deviation- nisp')


%%% Generating PDF at centre of domain
lo = 81;
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


for p = 1:n_outpce
U_coeff(p) = U(lo,p);    
U_pdf = U_pdf + U(lo,p).*psi_pdf(:,p);

end

[f_NI_5,xi_NI_5] = ksdensity(U_pdf);
load MCS_50000_03.mat;

figure(n_outpce+2)
plot(xi_NI_5,f_NI_5)
legend('Non-Intrusive -Quadrature')


% % mu_U_point = mean(U_pdf)
% % sigma_U_point = sqrt(var(U_pdf))
% % skew = skewness(U_pdf)
% % kur = kurtosis(U_pdf)

format long
U_coeff



