
function [kappa] = GetKLEterms_Normalized(x,d,t,file,s_spectral,ipce)

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
        % To get the same multipliers as specified in ajit's code , recheck
        % properly again
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


%        L_0 = exp(mu + 0.5 * sigma^2);

% Direct value to check with ajits code
        L_0 = 1.046;
%       L_0 = 0.04;


% Better to change loading everytime 
% filename = strcat('mindex',s_spectral.mindexname,'.dat');
m_index_in = s_spectral.mindex.m_index_out;

% % %  cd(cfolder) 
% filename = strcat('norm_squared',s_spectral.normname);
norm = s_spectral.norm;
    

% for npcin = 1:s_spectral.n_inpce
    
    numerator = 1;
    for ii = 1:dim
        numerator = numerator * gx_2D(ii,:)^(m_index_in(ipce,ii));
    end
        kappa = L_0*(numerator/sqrt(norm.norm_squared(ipce)));
        
        
     
% end 
  
  
 
  
   %..............Code from fortran modified .......................

% % %             xdash = x-a;
% % %                             
% % % %             4-Dimensional (L=4)
% % %             g(1) = sigma*multiplier_x(1)*multiplier_y(1)*(cos(omega(1)*xdash(1))*cos(omega(1)*xdash(2)));
% % %             g(2) = sigma*multiplier_x(1)*multiplier_y(2)*(cos(omega(1)*xdash(1))*sin(omega(2)*xdash(2)));
% % %             g(3) = sigma*multiplier_x(2)*multiplier_y(1)*(sin(omega(2)*xdash(1))*cos(omega(1)*xdash(2)));
% % %             g(4) = sigma*multiplier_x(2)*multiplier_y(2)*(sin(omega(2)*xdash(1))*sin(omega(2)*xdash(2)));
% % % % % % 
% % % % % %             switch ipce
% % % % % %                 case 1                           % 0th-Order (p=0)
% % % % % %                   y = 1;
% % % % % %                 case 2                           % 1st-Order (p=1)
% % % % % %                   y = g(1);
% % % % % %                 case 3
% % % % % %                   y = g(2);
% % % % % %                 case 4
% % % % % %                   y = g(3);
% % % % % %                 case 5
% % % % % %                   y = g(1)^2/2.0;
% % % % % %                 case 6                           %2nd-Order (p=2)
% % % % % %                   y = g(1)*g(2);  
% % % % % %                 case 7
% % % % % %                   y = g(1)*g(3);
% % % % % %                 case 8
% % % % % %                   y = g(2)^2/2; 
% % % % % %                 case 9
% % % % % %                   y = g(2)*g(3); 
% % % % % %                 case 10   
% % % % % %                   y = g(3)^2/2.0;   
% % % % % %       
% % % % % %             end 
% % % % % %        
% % % % % %             kappa = L_0*y;

%             gx_2D = 0;



 %%%%
% % %  fprintf('Mesh points are')
% % %   x(1)
% % %   x(2)
% % %   fprintf('Eigen function value using og method')
% % %   gx_2D
% % %   fprintf('Eigen function value using Ajits method')
% % %   g
% % %   
  


end