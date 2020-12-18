%% program to create a log-normal permeability field suitable for
%representing multiscale heterogeneity in porous media
%%

clc
clear variables

%    Defining the grid structure
x = 1:8:64;
y = 1:8:64;
[X,Y] = meshgrid(x,y);
[XX,YY] = meshgrid(1:32);

%    Initializing constants \beta 
beta1 = 0.5;
beta2 = 0.7;

n = size(x,2);

a = zeros(n^2,2);
for i=1:n
    for j=1:n
        a(n*(i-1)+j,:) = [X(i,j)  Y(i,j)];
    end
end
%    Defining the covariance function for Gaussian field
%c =1;
c     = 1/(1 - 2^(-beta1));
Sigma = c*ones(n*n);
for k =1:n*n
    for l=1:n*n
        if k~=l
             %Sigma(k,l) = (sum((a(k,:)-a(l,:)).^2))^(0.5);
             Sigma(k,l) = c*(sum((a(k,:)-a(l,:)).^2))^(-beta1/2);
        end
    end
end

                
%    Simulating Gaussian field with mean = 0 and covariance = sigma
R = mvnrnd(zeros(1, n*n),Sigma, n*n); R32 = R(1:32,1:32);
%    Using pseudocolor plotting function to plot the Gaussian field
figure(50); pcolor(R32); shading interp; colormap jet; colorbar; title('Gaussian field')
export_fig 'Kraster32new.pdf' '-opengl'
%    Computing permeability from the Gaussian field
k0 = 3.4253; varphi = 1;
KK = k0*exp(varphi*R); KK = 100*KK; KK = KK(1:32,1:32); 
save K32new.mat
% Plotting permeability
figure(51)
surf(XX,YY,log(KK)); shading interp; colorbar; title('Random permeability field')
xlim([0 32]); ylim([0 32]); colormap jet;
export_fig 'KK32new.pdf' '-opengl'