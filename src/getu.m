function out = getu(A,B)
maxit = 600; tic
out =bicgstab(A,B,9e-4,maxit); toc
% out=A\B;