

%%% Program to plot the standard deviation in PC coefficients from NISP using MC 
%%% with increasing sample size, both batch and quadrature calculated
clearvars
load PC_coeff_std.mat
load PC_coeff_quad.mat
load batch_sd.mat
load batch_mean.mat

lo = 9;
x = 10000:1000:50000;
y = 10000:10000:50000;

figure(1)
for n_pc = 1:8
    
    for i = 1:size(x,2)
        
        pc_coeff_COV = pc_coeff_std(lo,n_pc)/abs(pc_coeff(lo,n_pc));
        error_sd(n_pc,i) = pc_coeff_COV/sqrt(x(i));

    end
  
 % error_sd(n_pc,:) =  error_sd(n_pc,:)./max(error_sd(1,:));
    plot(x,error_sd(n_pc,:))
    hold on
end

for n_pc = 1:8
batch_cov(n_pc,:) = batch_sd(n_pc,:)/abs(pc_coeff(lo,n_pc)); 
% batch_cov(n_pc,:) = batch_sd(n_pc,:)./abs(batch_mean(n_pc,:)); 
plot(y,batch_cov(n_pc,:),'*')
hold on
end









