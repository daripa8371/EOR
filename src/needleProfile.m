%% Code to find needle lift profile
%% Rohit Mishra

% load('xanthane_prop.mat');
% f = fit( xanthane_prop(1:7,1), xanthane_prop(1:7,3),'rat23' );
% 
% 
% plot( f, xanthane_prop(:,1), xanthane_prop(:,3 ));

x = 0:1e-5:4e-3;

p1 =  0 ; %(-2.162e+23, -1.808e+23)
       p2 = 0;%  4.438e+21 ;% (4.029e+21, 4.846e+21)
       p3 =  0;%-4.352e+19 ;% (-4.766e+19, -3.938e+19)
       p4 =   5.54826989e13;%2.457e+17 ; %(2.216e+17, 2.699e+17)
       p5 =   - 7.04908939e11   ;% (-9.702e+14, -7.913e+14)
       p6 =   3.46202548e9 ; %(1.86e+12, 2.296e+12)
      p7 = - 8.29688306e6;
      p8 = 9871.79147;
      p9 = - 4.6569122;
      p10 = 2.90700300e-11;
          mult = p1*(x.^9) + p2*(x.^8) + p3*(x.^7) + p4*(x.^6) + ...
                    p5*(x.^5) + p6*(x.^4) + p7*(x.^3) + p8*(x.^2) + p9*x + p10;

plot(x,mult);
x1 = 0.0014632;
mult1 = p1*(x1^9) + p2*(x1^8) + p3*(x1^7) + p4*(x1^6) + ...
                    p5*(x1^5) + p6*(x1^4) + p7*(x1^3) + p8*(x1^2) + p9*x1 + p10;

% 
% Column1 = 0.00135:1e-7:0.004;
% Column2=2.90700300e-11 - 4.6569122*Column1 + 9871.79147*Column1.^2 ...
%     - 8.29688306e6*Column1.^3 + 3.46202548e9*Column1.^4 ...
% - 7.04908939e11*Column1.^5 + 5.54826989e13*Column1.^6;
% 
% 
% plot(Column1,Column2);
% hold on;
% plot(data(:,1), data(:,2 ));
%mult = p1*pow(t,9) + p2*pow(t,8) + p3*pow(t,7) + p4*pow(t,6) + p5*pow(t,5) + p6*pow(t,4) + p7*pow(t,3) + p8*pow(t,2) + p9*t + p10