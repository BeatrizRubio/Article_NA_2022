function BDA =BDAGram_matrix(n,m)

%Bidiagonal decomposition of the (n+1)x(n+1) Gram matrix of 
%Bernstein basis of negative degree m.

BDA=zeros(n+1);


%Computation of the multipliers m_{i,j}

for i=2:n+1
     for j=1:i
         BDA(i,j)=(m+i-2)*(2*m+i-3)/((2*m+i+j-3)*(2*m+i+j-4));
     end
end  

%Computation of the pivots p_{i,i}
 for i=1:n+1
     BDA(i,i)=(nchoosek(m+i-2,i-1)*factorial(i-1)*factorial(2*m+i-3))^2/(factorial(2*m+2*i-4)*factorial(2*m+2*i-3));
 end

%Computation of the multipliers tilde m_{i,i}
 
for i=1:n+1
    for j=i+1:n+1
        BDA(i,j)=(m+j-2)*(2*m+j-3)/((2*m+i+j-3)*(2*m+i+j-4));
    end
end


