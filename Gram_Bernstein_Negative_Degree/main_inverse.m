clear all
clear all
format longE
%Inverse matrix

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein of Negative Degree  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_Inverse.nb

n=24;
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

%Inverse Matrix

IB=TNInverseExpand(BDA);
IM=inv(A);
dlmwrite('inverseGramB.csv',IB,'precision','%.45f');
dlmwrite('inverseGramM.csv',IM,'precision','%.45f');

 %function A=TNInverseExpand(B)  
%Computes directly the inverse a square TN matrix whose bidiagonal
% bidiagonal decomposition is stored in B, using the 
% results on the factorization of A and its inverse presented in:
% Ana Marco, Jose-Javier Martinez:  Accurate computations with totally 
% positive Bernstein-Vandermonde matrices.
% Electronic Journal of Linear Algebra, Volume 26 (2013): 357--380.



%




