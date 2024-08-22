
function [kappa] = GetKLEterms(x,d,t,file,s_spectral,xi)

% Finding the eigen functions using tensor product in the order of
% descending lambda

%%% Code is actually being called to generate same values 9 times can be
%%% improved

omega = file.omega;
lambda_multi = file.lambda_multi;
lambda = file.lambda;
a = file.a;
b = file.b;


dim = s_spectral.dim;

mu = s_spectral.mu_g;
sigma = s_spectral.sigma_g;

% % % eigfun_x = zeros(length(omega),1);
% % % eigfun_y = zeros(length(omega),1);


for i = 1:dim
   
    num_1 = lambda_multi(i,2);
    num_2 = lambda_multi(i,3); 
    
    if(mod(num_1,2)~=0) % Odd
        den_odd = sqrt(a + sin(2*omega(num_1)*a)/(2*omega(num_1)));
        eigfun_x(i,:) = cos(omega(num_1).*(x(1)-a))/den_odd;
        % To get the same multipliers as specified in ajit's code
        multiplier_x(num_1) = sqrt(lambda(num_1))/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(num_1)*a)/(2*omega(num_1)));
        eigfun_x(i,:) = sin(omega(num_1).*(x(1)-a))/den_even;
        multiplier_x(num_1) = sqrt(lambda(num_1))/den_even;
    end
    
    
    if(mod(num_2,2)~=0) % Odd
        den_odd = sqrt(a + sin(2*omega(num_2)*a)/(2*omega(num_2)));
        eigfun_y(i,:) = cos(omega(num_2).*(x(2)-a))/den_odd;
        multiplier_y(num_2) = sqrt(lambda(num_2))/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(num_2)*a)/(2*omega(num_2)));
        eigfun_y(i,:) = sin(omega(num_2).*(x(2)-a))/den_even;
        multiplier_y(num_2) = sqrt(lambda(num_2))/den_even;
    end
    
    
    eigfun_xy(i,:) = eigfun_x(i,:).*eigfun_y(i,:);
    gx_2D(i,:) = sigma.*(sqrt(lambda_multi(i,1))).*eigfun_xy(i,:);
end


%   L_0 = exp(mu + 0.5 * sigma^2);

% Direct value to check with ajits code
    L_0 = 0.04;


% Better to change loading everytime 
% filename = strcat('mindex',s_spectral.mindexname,'.dat');
m_index_in = s_spectral.mindex;

% % %  cd(cfolder) 
% filename = strcat('norm_squared';s_spectral.normname);
norm = s_spectral.norm;
    

 for npcin = 1:s_spectral.n_inpce
    
    numerator = 1;
    for ii = 1:dim
        numerator = numerator * gx_2D(ii,:)^(m_index_in(npcin,ii));
    end
        kappa_coeff(npcin) = L_0*(numerator/norm.norm_squared(npcin));
           
 end 
  
  

%%% Multiplying chaos terms for MCS
 
switch (s_spectral.dim)

         case 2
                lterms(1) = 1.0;
                lterms(2) = xi(1);
                lterms(3) = xi(2);
                lterms(4) = (xi(1)^2 - 1);
                lterms(5) = xi(1)*xi(2);
                lterms(6) = (xi(2)^2 - 1);
                lterms(7) = (xi(1)^3 - 3*xi(1));
                lterms(8) = (xi(1)*xi(1)*xi(2) - xi(2));
                lterms(9) = (xi(2)*xi(2)*xi(1) - xi(1));
                lterms(10) = (xi(2)^3 - 3*xi(2));

         case 3
                lterms(1) = 1.0;
                lterms(2) = xi(1);
                lterms(3) = xi(2);
                lterms(4) = xi(3);
                lterms(5) = (xi(1)^2 - 1);
                lterms(6) = xi(1)*xi(2);
                lterms(7) = xi(1)*xi(3);
                lterms(8) = (xi(2)^2 - 1);
                lterms(9) = xi(2)*xi(3);
                lterms(10) = (xi(3)^2 - 1);
                lterms(11) = (xi(1)^3 - 3*xi(1));
                lterms(12) = (xi(1)^2*xi(2) - xi(2));
                lterms(13) = (xi(1)^2*xi(3) - xi(3));
                lterms(14) = (xi(2)^2*xi(1) - xi(1));
                lterms(15) = xi(1)*xi(2)*xi(3);
                lterms(16) = (xi(3)^2*xi(1) - xi(1));
                lterms(17) = (xi(2)^3 - 3*xi(2));
                lterms(18) = (xi(2)^2*xi(3) - xi(3));
                lterms(19) = (xi(3)^2*xi(2) - xi(2));
                lterms(20) = (xi(3)^3 - 3*xi(3));
end
 
kappa = 0;
  
for n_l = 1:s_spectral.n_inpce

kappa = kappa + kappa_coeff(n_l)*lterms(n_l);

end


   %..............Code from fortran modified .......................

% % %             xdash = x-a;
% % %                             
% % % %             4-Dimensional (L=4)
% % %             g(1) = sigma*multiplier_x(1)*multiplier_y(1)*(cos(omega(1)*xdash(1))*cos(omega(1)*xdash(2)));
% % %             g(2) = sigma*multiplier_x(1)*multiplier_y(2)*(cos(omega(1)*xdash(1))*sin(omega(2)*xdash(2)));
% % %             g(3) = sigma*multiplier_x(2)*multiplier_y(1)*(sin(omega(2)*xdash(1))*cos(omega(1)*xdash(2)));
% % %             g(4) = sigma*multiplier_x(2)*multiplier_y(2)*(sin(omega(2)*xdash(1))*sin(omega(2)*xdash(2)));

  


end