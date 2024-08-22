
clc
% Moments of PCE - Multiplication Tensor 
clearvars;
syms x 
syms one

dim = 3;

in_ord = 2;

out_ord = 3;


% One dimensional polynomial order (maximum of input and output order)
if(in_ord > out_ord)
    ord = in_ord;
else
    ord = out_ord;
end

% Number of terms in input expansion (k index for MultiDimensional Moment )

n_pce_in = factorial(in_ord + dim)/(factorial(in_ord)*factorial(dim));

% Number of terms in output expansion (i and j index for MultiDimensional Moment)

n_pce_out = factorial(out_ord + dim)/(factorial(out_ord)*factorial(dim));


% 1-st order Hermite PC: it is always fixed to psi(1) = 1
psi(1) = one;

if(ord > 1 || ord == 1 )
% 2nd order Hermite PC: it?s always fixed to psi(2) = x if (nord > 0)
psi(2) = x;

% 3rd and more order Hermite PC?s are solved by recursive formula for i = 3 : nord+1
    for i = 3:ord+1
        psi(i) = x * psi(i-1) - (i-2) * psi(i-2);
    end
end

% 1D Moments of Hermite PC basis <psi_i psi_j psi_k>
moment_1D = zeros(ord + 1, ord + 1, ord + 1);

for k = 1 : ord + 1
    for j = 1 : ord + 1
        for i = 1 : ord + 1
             
                  if(k ==1 || j == 1 || i == 1)
                      psi = subs(psi,one,1);
                  end
                  mult = (1/(sqrt(2*pi)))*int( psi(i) *  psi(j)* psi(k) * exp(-(x)^2/2),-inf,inf);
                  mult = double(mult);
                  moment_1D(i,j,k) = mult;     %   *((1/(factorial(i-1))^(1/dim))); 
             
        end
    end
end

moment_1D;

% Creating Multi indices from UQTk

% % % point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);

[status,out] = system(['./gen_mi -p', num2str(in_ord), ' -q', num2str(dim) , '-x TO'])

m_index_in = load('mindex.dat')


[status,out] = system(['./gen_mi -p', num2str(out_ord), ' -q', num2str(dim) , '-x TO'])

m_index_out = load('mindex.dat')

% Multi Dimensional Moments from 1D moments <Psi_i Psi_j Psi_k> = 
  
    
    for i = 1:n_pce_in
        
        for j = 1: n_pce_out
            
            for k = 1:n_pce_out
                
                product = 1;
                norm_product = 1;
                
                    for d = 1: dim
                        
                       mi1 = m_index_in(i,d)+1;
                       mi2 = m_index_out(j,d)+1;
                       mi3 = m_index_out(k,d)+1;
                       moment = product * moment_1D(mi1,mi2,mi3);
                       product = moment;
                       
%                        if(d < 2 || d ==2)
                           norm_mik = m_index_out(k,d);
                           norm_moment = norm_product * factorial(norm_mik);
                           norm_product = norm_moment;
%                        end

                    end
                    
                  norm_squared(k) = norm_moment;
                  moment_NonNormalized(i,j,k) = moment;
                  moment_ND(i,j,k) = moment/norm_squared(k);
                  
            end
        end
    end
   
   moment_NonNormalized
   moment_ND ;
   norm_squared;
   
cd ../../../ssfem_matlab/ssfem_LNP_Intrusive/;


str1 = strcat('0',int2str(ord));

if dim <= 10
    str2 = strcat('000',int2str(dim));
else
    str2 = strcat('00',int2str(dim));   
end

str3 = strcat('Cijk',str1,str2);
str4 = strcat('norm_squared',str1,str2);
        
save(str3,'moment_ND');
save(str4,'norm_squared');

%save('/Users/sudhipv/documents/UQTk_Matlab/norm_k.mat','norm_squared');


% Quadrature weights and points

% % point to UQTk directopry
% uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-build/bin';
% cd(uqtk_path);


% Non-dimentionalized, 1-dimensional Hermite polynomials 
% 1-dimensional Hermite polynomials at quad-points