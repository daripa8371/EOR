
function out = z_func(x,y)
% out=(x-0.5)^2+(y-0.5)^2-0.15;
% 
% [th,r] = cart2pol(x,y);
% if th > -pi && th <pi/7-pi/5
% 	th = th+2*pi;
% end
% tht = pi/5;
% thr = pi/7;
% R = 1;
% for ii = 1:5
% 	if th>=thr+pi*(2*ii-2)/5 && th<thr+pi*(2*ii-1)/5
% 		out = R*sin(tht/2)/sin(tht/2+th-thr-2*pi*(ii-1)/5)-r;
% 		return;
% 	end
% 	if th>=thr+pi*(2*ii-3)/5 && th<thr+pi*(2*ii-2)/5
% 		out = R*sin(tht/2)/sin(tht/2-th+thr+2*pi*(ii-1)/5)-r;
% 		return;
% 	end
% end

%  out=(x)^2+(y)^2-0.015;
out=-(-(-x+y)-3);