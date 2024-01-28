function BDA =BDAGram_matrix(n)

%Bidiagonal decomposition of the Gram matrix of 
%Bernstein basis 

BDA=zeros(n+1);


%Computation of the multipliers m_{i,j}

for i=2:n+1
     for j=1:i
         BDA(i,j)=(n-i+2)*(2*n-i+3)/((2*n-i-j+3)*(2*n-i-j+4));
     end
end  

%Computation of the pivots p_{i,i}
for i=1:n+1
    BDA(i,i)=(nchoosek(n,i-1)*factorial(i-1)/factorial(2*n-i+2))^2*factorial(2*n-2*i+2)*factorial(2*n-2*i+3);
end

%Computation of the multipliers tilde m_{i,i}
 
for i=1:n+1
    for j=i+1:n+1
        BDA(i,j)=(n-j+2)*(2*n-j+3)/((2*n-i-j+3)*(2*n-i-j+4));
    end
end


