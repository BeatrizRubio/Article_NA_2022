clear all
format longE
%Inverse matrix

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_Inverse.nb

n=24;


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








