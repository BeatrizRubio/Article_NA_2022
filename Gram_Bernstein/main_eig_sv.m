clear all
format longE

%Eigenvalues and Singular Values

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_EV_SV.nb

n=19;


A=zeros(n+1);
%Gram matrix of  Bernstein  basis  

for i=1:n+1
	for j=1:n+1
		A(i,j)=nchoosek(n,i-1)*nchoosek(n,j-1)*factorial(i+j-2)*factorial(2*n-i-j+2)/factorial(2*n+1); 
    end 
end


 
%Bidiagonal decomposition of Bersntein Gram matrix of 
%Geometric basis  

BDA=BDAGram_matrix(n);
 


%Eigenvalues
 
EVB=min(TNEigenValues(BDA))
EVM=min(eig(A))
dlmwrite('EVGramB.csv',EVB,'precision','%.45f');
dlmwrite('EVGramM.csv',EVM,'precision','%.45f');

%Singular values 
 
SVB=min(TNSingularValues(BDA));
SVM=min(svd(A));
dlmwrite('SVGramB.csv',SVB,'precision','%.45f');
dlmwrite('SVGramM.csv',SVM,'precision','%.45f');


% function a=TNEigenValues(B)
% Computes the eigenvalues of a TN matrix with bidiagonal decomposition
% stored in B
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written February 2003
%July 2015: Added the return of the eigenvector matrix V
%            Supported by SJSU's Woodward Fund

%function a=TNSingularValues(B);
%Computes the singular values of a TN matrix A with bidiagonal
% decomposition B=BD(A)
% Written February 2003
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.





