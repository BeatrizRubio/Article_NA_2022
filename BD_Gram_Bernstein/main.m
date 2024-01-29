clear all
format long E

%Bidiagonal decomposition  of Gram matrix of 
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

n=24;

A=zeros(n+1);

%Gram matrix of  Bernstein  basis  
for i=1:n+1
	for j=1:n+1
		A(i,j)=nchoosek(n,i-1)*nchoosek(n,j-1)*factorial(i+j-2)*factorial(2*n-i-j+2)/factorial(2*n+1); 
    end 
end

 
%Bidiagonal decomposition of Bersntein Gram matrix

 BDA=BDAGram_matrix(n)

%Linear system Ax=b
 b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29, -41,...
 46, -16, 74, -1, 2, -6, 7, -5, 1, -2, 6, -7, 5];
SolB=transpose(TNSolve(BDA,transpose(b)))
SolM=A\transpose(b)

%Inverse Matrix

IB=TNInverseExpand(BDA)

%Eigenvalues
 
EVB=min(TNEigenValues(BDA))


%Singular values 
 
SVB=min(TNSingularValues(BDA))




