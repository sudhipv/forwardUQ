

%% Exponential Covariance - Two dimension

% ......Sudhi Sharma P V ..............

%%Taking the tensor Product in Both directions


clear all
% delete 'eigfun_2D.dat'
b = 1;
a = 0.5; % Half the domain length 
mult = 0.5/a;


% Domain

square

x = p(1,:);
y = p(2,:);

KLE_dim = 3; 
mu_g  = 0;
sigma_g = 1;
ord_in = 2;


% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + KLE_dim)/(factorial(ord_in)*factorial(KLE_dim));


n = 1000;
wroot = zeros(1,n);
init = 1;
j_o = 1;
check_omega = 0;
step = 0;

% Odd terms - Equations are swapped in Roger Ghanem Book
for i= 1:n
    
fun =@(w)(1/b - w*tan(w*a));

[wroot(i) fval(i) flag] = fzero(fun,init);

 
        if(flag ~=1)
            init = init +0.01;
            continue
        end

       if(any(abs(check_omega - wroot(i))>0.01,1))
        
           omega_o(j_o) = wroot(i);
           check_omega(j_o) = omega_o(j_o);
           lambda_o(j_o) = 2*b/(1+b^2*omega_o(j_o)^2);
           init = omega_o(j_o);
           
           if(j_o>1)
                  step = omega_o(j_o) - omega_o(j_o-1);
           else
                  step = 2*pi*mult;
           end
       
           j_o = j_o+1;
           init = init +step;
        
       end
       
    
end

lambda_o
omega_o


% Even Terms
init = 3.60;
check_omega = 0;
j_e = 1;
step =0;
ne= 1000;


for i= 1:ne
    
fun_e =@(w)(w + (1/b)*tan(w*a));


[wroot_e(i) fval_e(i) flag_e] = fzero(fun_e,init);
 
        
% Condition added since froot gives a flag of 1 even after giving singular
% root value
if(wroot_e(i)>1)
    
    if(flag_e ~=1)
            init = init +0.1;
            continue
    end
        
         if(any(abs(check_omega - wroot_e(i))>0.01,1))
             omega_e(j_e) = wroot_e(i);
             check_omega(j_e) = omega_e(j_e);
             lambda_e(j_e) = 2*b/(1+b^2*omega_e(j_e)^2);
             init = omega_e(j_e);
             step = 2*pi*mult;
             j_e = j_e+1;
             init = init + step;
         else
            init = init +0.1; 
         end
         
else
    init = init +0.1;
end

    
end

lambda_e  
omega_e

% Ordered Omega and Lambda according to odd and even 
n = 10;
n_e = 1;
n_o = 1;
del_odd = 0;
del_eve = 0;

for i= 1:n
       
   if(mod(i,2) ~=0)
       lambda(i) = lambda_o(n_o);
       omega(i) = omega_o(n_o);
       n_o = n_o + 1;
   
   else
       lambda(i) = lambda_e(n_e);
       omega(i) = omega_e(n_e);
       n_e = n_e + 1;

   end
       
end
    
 lambda  = sort(unique(lambda),'descend')
 omega   = sort(unique(omega))


figure(3)
plot(1:length(lambda),lambda,'*')
ax = gca;
set(gca, 'YScale', 'linear','ytick', [0.0001 0.001 0.01 0.1 1]);


k = 1;
for i = 1:length(lambda)
    for j= 1:length(lambda)
        
        lambda_n(k,1) = lambda(i)*lambda(j);
        lambda_n(k,2) = i;
        lambda_n(k,3) = j;
        k = k+1;
    end 
end
[B,I] = sort(lambda_n(:,1),'descend');
lambda_multi = lambda_n(I,:,:)

% Finding the eigen functions using tensor product in the order of
% descending lambda

% % x = -0.5:0.01:0.5;
% % y=-0.5:0.01:0.5;
 eigfun_xy = zeros(KLE_dim,length(p));
 


for i = 1:KLE_dim
   
    num_1 = lambda_multi(i,2);
    num_2 = lambda_multi(i,3); 
    
    if(mod(num_1,2)~=0) % Odd
        den_odd = sqrt(a + sin(2*omega(num_1)*a)/(2*omega(num_1)));
        eigfun_x(i,:) = cos(omega(num_1).*(x-a))/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(num_1)*a)/(2*omega(num_1)));
        eigfun_x(i,:) = sin(omega(num_1).*(x-a))/den_even;
    end
    
    
    if(mod(num_2,2)~=0) % Odd
        den_odd = sqrt(a + sin(2*omega(num_2)*a)/(2*omega(num_2)));
        eigfun_y(i,:) = cos(omega(num_2).*(y-a))/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(num_2)*a)/(2*omega(num_2)));
        eigfun_y(i,:) = sin(omega(num_2).*(y-a))/den_even;
    end
    
    
    eigfun_xy(i,:) = eigfun_x(i,:).*eigfun_y(i,:);
    g_x_2D(i,:) = (sqrt(lambda_multi(i,1))*sigma_g^2).*eigfun_xy(i,:);
end



file = 'eigfun_2D.dat';
dlmwrite(file,g_x_2D);
% Number of terms in Expansion
% % L = 10; % Number of terms in PCE +1 : with 2nd order and 3 RV 
% +1 is because we are starting from 1 here

% % % file = 'eigfun_2D.dat';
% % % for i = 1:L
% % % eigfun_multi(:,:,i) = eigfun(lambda_multi(i,2),:)'*eigfun(lambda_multi(i,3),:);
% % % g_x_2D(:,:,i) = sqrt(lambda_multi(i,1))*eigfun(lambda_multi(i,2),:)'*eigfun(lambda_multi(i,3),:);
% % % dlmwrite(file,g_x_2D,'-append')
% % % g_x_2D(:,:,i) = g_x_2D(:,:,i).*sigma_g^2;
% % % end
 


for i=1:KLE_dim
figure(i)
pdesurf(p,t,eigfun_xy(i,:)')
shading faceted
xlabel('X')
ylabel('Y')
zlabel('\surd(\lambda_i_x * \lambda_j_y) g_i(x) g_j(y)')
end



%%% Lognormal PCE coefficients 2Dimension

% L = 10;
% % % g_0 = ones(length(x),length(y)).*mu_g;

% % % % g_i = zeros(length(x),length(y));
% % % % 
% % % % for i = 1:KLE_dim
% % % %     
% % % %     g_i = g_i + 0.5 * g_x_2D(:,:,i).^2;
% % % %     
% % % %     
% % % % end
% % % % 
% % % % L_o = exp(g_0 + g_i);

L_0 = exp(mu_g + 0.5 * sigma_g^2);


% PCE coefficients of Lognormal Process with inherent gaussian process


% % % L_x_2D(1,:) = L_0 .* ones(1,length(x)); % 1
% % % L_x_2D(2,:) = L_0 .* g_x_2D(1,:); % xi1
% % % L_x_2D(3,:) = L_0 .* g_x_2D(2,:); % xi2
% % % L_x_2D(4,:) = L_0 .* g_x_2D(3,:); % xi3
% % % L_x_2D(5,:) = L_0 .* g_x_2D(1,:).^2/2; % xi1^2-1
% % % L_x_2D(6,:) = L_0 .* g_x_2D(1,:).*g_x_2D(2,:); % xi1*xi2;
% % % L_x_2D(7,:) = L_0 .* g_x_2D(1,:).*g_x_2D(3,:); % xi1*xi3; 
% % % L_x_2D(8,:) = L_0 .* g_x_2D(2,:).^2/2; % xi2^2-1; 
% % % L_x_2D(9,:) = L_0 .* g_x_2D(2,:).*g_x_2D(3,:); % xi2*xi3; 
% % % L_x_2D(10,:) = L_0 .* g_x_2D(3,:).^2/2; % xi3^2-1; 
% % % 


m_index_in = load('mindex030003.dat');
norm = load('norm_squared040003.mat');
    

for npcin = 1:n_inpce
    
    numerator = 1;
    for ii = 1:KLE_dim
        numerator = numerator .* g_x_2D(ii,:).^(m_index_in(npcin,ii));
    end
        L_x_2D(npcin,:) = L_0*(numerator/norm.norm_squared(npcin));
        
        
     
end 

for i = 1:n_inpce
figure(10+i)
pdesurf(p,t,L_x_2D(i,:)')
shading faceted
end
legend



