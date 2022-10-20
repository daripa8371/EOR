function out=minmod(a,b)


out = ((abs(a)<=abs(b)).*(a.*b>0)).*a+((abs(b)<abs(a)).*(a.*b>0)).*b;
