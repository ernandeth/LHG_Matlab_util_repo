function U = makeUnfoldMat(S, psi);
% function U = makeUnfoldMat(S, psi);
% 
% create the "unfolding" matrix for SENSE reconstruction by following resipe
% from Pruessmann papers
%
%  U = pinv(SH * psi_inv *S)  * SH *psi_inv;
% 
% S is a matrix with sens. maps for all coils
% psi is the covariance matrix among the coils
%

Ncoils= size(psi,1);

psi_inv = pinv(psi);
SH = transpose(conj(S));
% SH = S';
U = pinv(SH * psi_inv *S)  * SH *psi_inv;

return