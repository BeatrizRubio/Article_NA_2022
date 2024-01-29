clear all
format longE

%Bidiagonal decomposition  of Gram matrix of 
%Bernstein of Negative Degree  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 


n=24;
m=10;

A=zeros(n+1);
%Gram matrix of  Bernstein  basis of Negative Degree m

for i=1:n+1
	for j=1:n+1
		A(i,j)=(nchoosek(m+i-2,i-1)*nchoosek(m+j-2,j-1)*factorial(i+j-2)*factorial(2*m-2))/factorial(2*m+i+j-3); 
    end 
end


%Bidiagonal decomposition of Gram matrix of Bernstein basis of Negative 
%Degree


BDA=BDAGram_matrix(n,m);

%Linear system Ax=b
 b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29, -41,...
 46, -16, 74, -1, 2, -6, 7, -5, 1, -2, 6, -7, 5];
SolB=transpose(TNSolve(BDA,transpose(b)))
SolM=A\transpose(b)

%Inverse Matrix

IB=TNInverseExpand(BDA)
IM=inv(A)


%Eigenvalues
 
EVB=min(TNEigenValues(BDA))
EVM=min(eig(A))


%Singular values 
 
SVB=min(TNSingularValues(BDA))
SVM=min(svd(A))



