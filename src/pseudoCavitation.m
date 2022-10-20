limitedAlpha1 = 1;
Cv_ = 1;
xSol = 0.0015;
pdimless = 580:1e5;
xSol_p1 = exp(-6.0895+1.428*log(3) + 0.9759*log(pdimless/1e6));
xSol_p2 = exp(-4.9012+0.62131*log(3) + 0.98022*log(pdimless/1e6) - 0.0046896 *(pdimless/1e6));

mvCoeff1 = (limitedAlpha1)*Cv_*  ...
    (max(0,xSol - xSol_p1));
nDmax = log(3.5e14);
nDmin = log(2.3e11);
lognD = (-(nDmax - nDmin)/(xSol)*xSol_p2) + nDmax;
nD = exp(lognD);


% Cvmin = 0.1;
% Cvmax = 1;
% Cv =  (-(Cvmax - Cvmin)/(xSol)*xSol_p) + Cvmax;
pAmb = 1e5;
xSol_p1_amb = exp(-6.0895+1.428*log(3) + 0.9759*log(pAmb/1e6));
xSol_p2_amb = exp(-4.9012+0.62131*log(3) + 0.98022*log(pAmb/1e6)...
    - 0.0046896 *(pAmb/1e6));

% 
% plot (pdimless,xSol_p2,'--','LineWidth',1.5);
% figure(2);
% semilogy(pdimless,nD,'--','LineWidth',1.5);
% legend('Number density');


plot(pdimless,xSol_p1,'r--')
hold on;
plot(pdimless,xSol_p2,'b--')
% legend('N_2 solubility')