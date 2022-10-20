function out = getu(A,B)
maxit = 300; tic
out =bicgstab(A,B,[],maxit); toc
% out=A\B;