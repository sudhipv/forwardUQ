

function [] = GetLambdaSymmetric_2D(b,a,dim)

delete lambda_2D.mat

% b = 1;
% a = 0.5; % Half the domain length 
n = 50;
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
                  step = 2*pi;
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
ne= 1100;


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
             step = 2*pi;
             j_e = j_e+1;
             init = init + step;
         else
            init = init +0.1; 
         end
         
else
    init = init +0.1;
end

    
end

lambda_e;  
omega_e;

% Ordered Omega and Lambda according to odd and even 
n = 20;
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
    
 lambda  = sort(unique(lambda),'descend');
 omega   = sort(unique(omega))



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
lambda_multi = lambda_n(I,:,:);

figure(3)
plot(1:length(lambda_multi),lambda_multi(:,1),'*')
ax = gca;
set(gca, 'YScale', 'log');


figure(4)
total_sum = sum(lambda_multi,1);
previous = 0;
partial_lambda = zeros(1,length(lambda_multi));
for ss = 1:length(lambda_multi)
partial_lambda(ss) = (lambda_multi(ss,1)+previous)/total_sum(1);
previous = previous + lambda_multi(ss,1);
end
plot(1:length(lambda_multi),partial_lambda,'*')
ax = gca;
set(gca, 'YScale', 'log');
 
save('lambda_2D.mat','lambda_multi','omega','a','b','lambda')

end