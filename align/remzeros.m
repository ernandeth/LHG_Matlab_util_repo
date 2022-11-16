function out=remzeros(matrix,col_no)
%REMZEROS Used to remove zero rows from a matrix.
%Takes as input the matrix and the column in which to look 
%for zeros.  If a zero is found the row is removed.
%
%  Example   Y=remzeros(matrix,col_no)

matrix(matrix(:,col_no)==0,:)=[];
out=matrix;
return;
