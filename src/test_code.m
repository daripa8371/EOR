% load('SIR_20.mat');
% C1 = COC';
% load('SIR_30.mat');
% C2 = COC';
% load('SIR_40.mat');
% C3 = COC';
% load('DIR_1.mat');
% C4 = COC';
% load('DIR_2.mat');
% C5 = COC';
% 
% x = 1:12500;
% 
% plot(x,C1(:,1));
% hold on;
% plot(x,C2(:,1));
% hold on;
% plot(x,C3(:,1));
% hold on;
% plot(x,C4(:,1));
% hold on;
% plot(x,C5(:,1));
% hold on;
% 













% clear;
% %% Test
% mu = 0.05;
% sigma = .4;
% xx = zeros;
% yy = zeros;
% 
% for i=1:100
%     x = i/1000;
% 
% 
% pd = makedist('Normal',mu,sigma);
% % pd = makedist('Weibull','a',0.5,'b',1.8)
% y = pdf(pd,x)
% xx(i) = x;
% yy(i) = y;
% end
% 
% 
% plot(xx,yy);

% 
% wppm = 50:1:1750;
% epsilon = 0.0001763 *wppm.^2  -0.0194*wppm + 5.149;
% plot(wppm,epsilon);
% gamma_dot = 0.01;
% for i=1:5
%    
% wppm = 1:10:1750;
% epsilon = 0.0001763 *wppm.^2  -0.0194*wppm + 5.149;
% power_n = 1.928e-07*wppm.^2 + -0.0006604*wppm + 0.8995;
% miua = 1e-3*epsilon .* (gamma_dot.^(power_n-1));
% plot (wppm,miua);
% hold on;
% gamma_dot = gamma_dot+0.05;
% end


% c = 0:0.0001:0.1;
% rho_water = 1000;
%         rho_xanthane = 1500;
%         w1=rho_xanthane * c;
%         w2=rho_water*(1-c);
%         wppm = ((w1)./ (w1 +w2))*1e4;
% 
% plot(c,wppm);
% diff = miua_1_20 - miua_0_20;
% 
% contourf (xCopy,yCopy,diff,5);

% load('miua_mixed.mat');
% plot(miua(15,:));
% load('miua_unmixed.mat');
% hold on;
% plot(miua(15,:));
% legend('Mixed','Unmixed');
% 
% load('procSpeed');
% plot(data1(:,1),data1(:,2),'--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
% 
% load('axial_voidInBeam.mat');
% 
% plot(exp_axial(:,1),exp_axial(:,2),'*');
% hold on;
% plot(num_axial(:,1),num_axial(:,2),'r--*');
%  hold on;
%  plot(lit_axial(:,1),lit_axial(:,2),'b--o');
%  legend('Experiment','Numerical','Battistoni et al.');
% 
% xlabel('x/L [-]'); 
% ylabel('Void in beam path [\mum]');

% load('pressureDrop.mat');
% 
% plot(pressure_exp(:,1),pressure_exp(:,2),'b--o');
% hold on;
% plot(pressure_num(:,1),pressure_num(:,2),'r--*' );
% legend('Experiment','Numerical');
% 
% figure(1);
% load('radial1_voidInBeam.mat');
% plot(exp_axial(:,1),exp_axial(:,2),'*');
% hold on;
% plot(num_axial(:,1),num_axial(:,3),'r--*');
%  hold on;
%  plot(lit_axial(1:27,1),lit_axial(1:27,2),'b--o');
%  legend('Experiment','Numerical','Converge-Battistoni');
% xlabel('r/R [-]'); 
% ylabel('Void in beam path [\mum]');
% 
% figure(2)
% load('radial7_voidInBeam.mat');
% 
% plot(exp_axial(:,1),exp_axial(:,2),'*');
% hold on;
% plot(num_axial(1:20,1),num_axial(1:20,3),'r--*');
%  hold on;
%  plot(lit_axial(1:22,1),lit_axial(1:22,2),'b--o');
%  legend('Experiment','Numerical','Converge-Battistoni');
% xlabel('r/R [-]'); 
% ylabel('Void in beam path [\mum]');


% load('stdaxial_voidInBeam.mat');
% 
% plot(exp_axial(:,1),exp_axial(:,2),'*');
% hold on;
% plot(num_axial(:,1),num_axial(:,3),'r--*');
% hold on;
% plot(lit_axial(:,1),lit_axial(:,2),'b--o');
% legend('Experiment','Numerical','Battistoni et al.');
% xlabel('x/L [-]'); 
% ylabel('Void in beam path [\mum]');
% 
% figure(2)
% load('stdradial7_voidInBeam.mat');
% 
% plot(exp_axial(1:27,1),exp_axial(1:27,2),'*');
% hold on;
% plot(num_axial(:,1),num_axial(:,2),'r--*');
%  hold on;
%  plot(lit_axial(1:22,1),lit_axial(1:22,2),'b--o');
%  legend('Experiment','Numerical','Converge-Battistoni');

% load('total_voidFraction.mat');
% 
% plot(exp_total(1:10,1),exp_total(1:10,2),'*');
% hold on;
% plot(exp_total(1:10,1),num_total(1:10,1),'r--*','LineWidth',2);
%  hold on;
%  plot(exp_total(1:10,1),num_total(1:10,2),'b--o','LineWidth',2);
%  legend('Experiment','Mixture Model','Multi-Fluid Model');

% 
% grid study
% load('gridStudy.mat');
% 
% plot(grid1(:,1),grid1(:,2),'b--o');
% hold on;
% plot(grid2(:,1),grid2(:,2),'r--*');
% hold on;
% plot(grid3(:,1),grid3(:,2),'-');
% hold on;
% plot(grid4(:,1),grid4(:,2),'*');
% legend('9\mum','18\mum','24\mum','Experiment');


% t = 0:1e-4:4e-3;
%         p4 =   5.54826989e13; % (4.307e+16, 4.485e+16)
%         p5 =    - 7.04908939e11; % (-1.33e+14, -1.274e+14)
%         p6 =   3.46202548e9; % (2.369e+11, 2.479e+11)
%         p7 =   - 8.29688306e6; % (-2.839e+08, -2.707e+08)
%         p8 =   9871.79147; % (1.74e+05, 1.831e+05)
%         p9 =      - 4.6569122; % (-51.03, -48.38)
%         p10 =   2.90700300e-11; % (-8.068e-09, 8.068e-09)
% 
% 
%         mult1 = p4*t.^6 + p5*t.^5 + p6*t.^4 ...
%                             + p7*t.^3 + p8*t.^2 + p9*t + p10;
%         
%                         
% plot(t,mult1);


% T = 600:10:1200;
% R = 1.9872;
% 
% P = R*T*2*8.3*10^13.*exp(-14413./(R*T))./(2.8*(10^18).*T.^(-0.86));
% 
% logP = log10(P*41.86);
% plot(T,logP,'*');
clear;
% src = 1;
% src_total = 0;
% rate =0.1;
% tcal = 0;
% src_save = zeros(101);
% while tcal<=100
%    src = src + rate*tcal*0.22941;
%     
%     
%     src_total = src_total + src;
%     src_save(tcal+1) = src;
%     tcal = tcal+1;
% end
% plot(src_save);
% clear;
% xx = -100:100
% xx0 = 1;
% v = (xx-xx0).^2;
% plot(xx,v)

% P = 0:100:10e6;
% x1=exp(-6.0895+1.4280*log(3)+0.9759*log(P/1e6));
% 
% plot (P,x1,'--','LineWidth',2);


%  a1=0.03298677e+02;
%         a2=0.14082404e-02;
%         a3=-0.03963222e-04;
%         a4=0.05641515e-07;
%         a5=-0.02444854e-10;
%         T = 300;
% cp = 8314/28*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

% k = 0.431638;
% Am = 0.1627;
% Tc = 126.192;
% T = 298;
% rho=398;
% MW=28.0134e-3;   
% Vm = MW/rho;
% Bm = 2.4e-5;
% G = k*sqrt(T/Tc)/(1+k*(1-sqrt(T/Tc)));
%      dAmdT = -(1/T)*Am*G;
%      R_molar = 8.3145;
%      Pc=3.3958e+06;
% d2Am_dT2=(0.457236/2)*(R_molar^2/Pc)*(k*(k+1))*(T/Tc).^(-3/2);     
%   dpdTv = (R_molar/(Vm-Bm))-(dAmdT/(Vm^2+2*Vm*Bm-Bm^2));
%      dpdvT = (-R_molar*T/(Vm-Bm)^2)*(1-2*Am*(R_molar*T*(Vm+Bm)*((Vm/(Vm-Bm))+(Bm/(Vm+Bm)))^2)^(-1));    


% 
%  load('miua_mixed.mat');
%  miuaMixed = miuaTcal;
%  load('miua_unmixed.mat');
%  miuaUnmixed = miuaTcal;
%  miuaDiff = miuaMixed - miuaUnmixed;
%  for i=1:6
%      figure(i);
%  contourf(miuaDiff(:,:,2),5);
%  end
%  


load('MFR.mat');
figure(1);
 plot(MFR(1:92,1),MFR(1:92,2)/112.48,'LineWidth',2);
hold on;
plot(MFR(1:94,3),MFR(1:94,4)/112.48,'r-o','LineWidth',2);
figure(2);
plot(MFR(:,5)*1e3,MFR(:,6)/44.19e-2,'LineWidth',2);

% DIFF = zeros;
% T=300:100:1600;
%  diff0 = 2.5563e-5;
%  Tref  = 273.11;
% Sref  = 110.56;
% for i=1:size(T,2)
%     DIFF(i) = diff0 * (T(i)/Tref).^1.5 * (Tref+Sref)/(T(i)+Sref);
% end
% 
% plot(T,DIFF);

