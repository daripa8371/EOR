%% Determine which subdomain the grid points belong to
%
% Evaluate the level set function on the grid points
% Negative output signifies the domain $$\Omega^-$$.
% Positive output signifies the domain $$\Omega^+$$.

function out=get_phi_test(para)


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
[jj,ii] = meshgrid(1:n+1,1:m+1);
get_phi = z_func_test(left+(ii-1)*dx, bottom+(jj-1)*dy);

%-------

out=get_phi;