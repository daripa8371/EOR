%% Code for post processing of finger width for QFS geometry(Rohit Mishra; 04/06/2022)
function [interfaceLoc,avgFW] = postProcess_fingerWidthQFS(UU)

count=0;
r=1;
x=1;
conc = zeros(30,1);
while(r<=30)
    while(x<=30)
        if(x>r)
            break;
        end
           
        y = fix(abs(sqrt(r^2 - x^2)));
        y = max(y,1);
        conc(r) = conc(r) + UU(x,y);
        count = count + 1;
        x=x+1;
    end
    conc(r) = conc(r)/count;
    count=0;
    x=1;
    r=r+1;
end
count=0;
conc=conc';
locMean = 0;
for i=1:30
  if (conc(i)>0.211 && conc(i)<0.99)
            locMean =  locMean + i;
           count = count+1;
  end
end
locMean = locMean/count;

interfaceLoc = fix(locMean);




        
r = interfaceLoc;
x = 1;
%% Find presence of saturation along the mixing layer
while(x<=30)
        if(x>r)
            break;
        end
           
        y = fix(abs(sqrt(r^2 - x^2)));
        y = max(y,1);
        if(UU(x,y) > 0.5)
            checkConc(x) = 1;
        else
            checkConc(x) = 0;
        end
        x=x+1;
end
x=1;
check = 0;
count = 1;
i = 1;
FW=0;
FWTotal = zeros;
countF = zeros;
checkConc = checkConc';
while(x<=size(checkConc,1))
    if(checkConc(x) == check)
        FW = FW+1;
        count = count + 1;
    else
        FWTotal(i) = FW;
        countF(i) = count;
        count = 1;
        i = i+1;
    end
    x=x+1;
end
FWTotal(i) = FW;
countF(i) = count;
avgFW = mean(FWTotal./countF);


end


