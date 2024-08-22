

%%% NON INTRUSIVE SSFEM
clearvars
%rng('default')
% Load the mesh
% square
% square_refined

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';

% Assemble matrix
A = AssembleMatrix(p, e, t, 'Poisson', [], 0);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson', [], 0);

% Solve the linear system
% U = A \ b;

n = 25000;
xi = randn(1,n);


psi(:,1) = ones(1,n);
psi(:,2) = xi;
psi(:,3) = xi.^2-1;
psi(:,4) = xi.^3-3.*xi;

mu_g = 0;
sigma_g = 0.3;

mu_l = exp(mu_g + sigma_g^2/2);
kappa(1) = mu_l;
kappa(2) = mu_l*sigma_g;
kappa(3) = mu_l*sigma_g^2/2;   
kappa(4) = mu_l*sigma_g^3/6;   

U = zeros(size(p,2),n);

for i = 1:n

   A_MCS = (A.*kappa(1).*psi(i,1)+A.*kappa(2).*psi(i,2)+A.*kappa(3).*psi(i,3));
 
%   A_MCS = A.*(exp(mu_g+sigma_g*xi(1,i))); +A.*kappa(4).*psi(i,4)
   
%%% Gaussian R V
%    A_MCS = A.*((mu_g+sigma_g*xi(1,i)));
   U(:,i) = A_MCS\b;     
       
end
  
% Mean and Standard Deviation of solution
 U_mean = sum(U,2)/n;
 
 sum = zeros(size(p,2),1);
 for j = 1:n
   
     sum = sum + (U(:,j) - U_mean).^2;
     
     
 end
 U_std = sum/(n-1);
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

%% 

loc = 81;
[f,xi] = ksdensity(U(loc,:));
% load MCS_50000_03.mat
figure(100)
plot(xi,f)
title('PDF of solution at centre of Domain KDE')
legend('MCS 50000','square- 289 nodes')


mu_U_point = mean(U(loc,:))
sigma_U_point = sqrt(var(U(loc,:)))
skew = skewness(U(loc,:))
kur = kurtosis(U(loc,:))

%%% Use below code to save the values for PDF generation saved in results
%%% folder.

% cd ../results;
% save('ksd_ssfem_MCS_10000.mat','xi_MCS','f_MCS')

