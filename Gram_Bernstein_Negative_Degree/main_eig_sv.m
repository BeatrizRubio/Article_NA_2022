clear all
format longE
%Eigenvalues and Singular values

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein of Negative Degree  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_EV_SV.nb

n=9;
m=10


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


%Eigenvalues
 
EVB=min(TNEigenValues(BDA));
EVM=min(eig(A));
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





