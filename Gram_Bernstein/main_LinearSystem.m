clear all
format rat
%Linear system Ax=b

%Experiments results of the
%Bidiagonal decomposition  of Gram matrix of 
%Bernstein  basis  Mass Matrix 
%E. Mainar, J.M. Peña, B. Rubio, 

%See experimental results in Mathematica: Gram_LinearSystem.nb

n=2;


A=zeros(n+1);
%Gram matrix of  Bernstein  basis  

for i=1:n+1
	for j=1:n+1
		A(i,j)=nchoosek(n,i-1)*nchoosek(n,j-1)*factorial(i+j-2)*factorial(2*n-i-j+2)/factorial(2*n+1); 
    end 
end
A

%  
% %Bidiagonal decomposition of Bersntein Gram matrix of 
% %Geometric basis  
% 
 BDA=BDAGram_matrix(n)
% 
 BDA2=TNBD(A)
%  
% 
% %Linear system Ax=b
% 
% %  b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29, -41,...
% %  46, -16, 74, -1, 2, -6, 7, -5, 1, -2, 6, -7, 5];
% % 
% % SolB=transpose(TNSolve(BDA,transpose(b)))
% % SolB2=transpose(TNSolve(BDA2,transpose(b)))
% %SolM=A\transpose(b)
% % 
%  %dlmwrite('sistemaGramB.csv',SolB,'precision','%.45f');
%  %dlmwrite('sistemaGramM.csv',SolM,'precision','%.45f');
% % 
% % %function TNSolve(B,b)
% % %Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)
% % 
% % %Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
% % %functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.
% % 
% % 
% % 
% % 
% % 
