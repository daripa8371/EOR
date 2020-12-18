%% Determine which subdomain the grid points belong to
%
% Evaluate the level set function on the grid points
% Negative output signifies the domain $$\Omega^-$$.
% Positive output signifies the domain $$\Omega^+$$.
% flag = 1 means Quarter five spot flood
% flag = 2 means HS flood
% flag = 3 means Five-spot parallel flood
% flag = 4 means Five-spot diagonal flood
function out=get_phi_test(para,flag)



m = para.box.m;
n = para.box.n;

dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;

% ------Loop Implementation
% get_phi=zeros(n+1,m+1);
% for ii=1:m+1
%     for jj=1:n+1
% 
%             get_phi(jj,ii)=z_func_test(left+(ii-1)*dx,bottom+(jj-1)*dy);
% 
%     end
% end
%-------------------

% --- Vectorized implementation
[ii,jj] = meshgrid(1:m+1,1:n+1);
if flag == 1
    get_phi = z_func_test(left+(ii-1)*dx, bottom+(jj-1)*dy);
elseif flag == 2
    get_phi = z_func_test_hs(left+(ii-1)*dx, bottom+(jj-1)*dy);
elseif flag == 3
    get_phi = z_func_test_parallel(left+(ii-1)*dx, bottom+(jj-1)*dy);
elseif flag == 4
    get_phi = z_func_test_diagonal(left+(ii-1)*dx, bottom+(jj-1)*dy);

end
%-------

out=get_phi;