%% Code for post processing of finger width (Rohit Mishra; 06/23/2019)
function [interface, meanFingerWidth,iterX] = postProcess_fingerWidth(UU)

interface = zeros (29,1);
%% add the mixing layer code to find the actual x position to find the width

% iterX = 1;
%              while (iterX <= 30)
%                  
%                  meanUU = mean(UU(:,iterX));
%                  if (meanUU<0.5)
%                      %iterX_save(tcal) = iterX;
%                      break
%                  end
%              iterX = iterX + 1;    
%                  
%              end
             
iterY = 1;
iterX = 1;
meanUUsave = zeros(29);
 countIter = 0;
    storeUU = 0;
%% Find the average concentration for each Y level
while (iterY <=29)
   
    while (iterX <=29)
        if (UU(iterY,iterX)>0.21 && UU(iterY,iterX)<0.99)
            countIter = countIter + 1;
            storeUU = UU(iterY,iterX) + storeUU;
            
        end
        iterX = iterX + 1;
    end
    meanUUsave(iterY) = storeUU/countIter;
     countIter = 0;
    storeUU = 0;
    iterY = iterY + 1;
    iterX=1;
end
iterY = 1;
iterX = 1;
while (iterY <=  29)
    while (iterX <=29)
        if (UU(iterY,iterX) < meanUUsave(iterY))
            
        interface(iterY,1) = iterX;
        break;
        end
        iterX = iterX+1;
    end
    iterX = 1;
    iterY = iterY + 1;
end


%% Find location of interface front

iterX = 1;
iterY = 1;
iterXSave = zeros(29);

while (iterX <= 29)
                 
 %meanUU = mean(UU(:,iterX));

 if (mean(UU(:,iterX))<mean(meanUUsave(:,1)))
     
     break
 end
iterX = iterX + 1;    

end




checkConc = zeros;

%% Find presence of saturation along the mixing layer

last = 0;
counter = 1;
meanFingerWidth = 0;
for ii=1:size(UU,1)
    for jj=1:size(UU,2)
        if UU(ii,jj)<meanUUsave(1,iterY)
            checkConc(ii,jj) = 1;
        else
            checkConc(ii,jj) = 0;
        end
    end
    
end


    totalConc = 0;
    jj = iterX-1;
  for i=1:1
      
        for ii=1:size(UU,2)
           % if jj == iterX-1 %|| jj == iterX || jj == iterX+1
                if checkConc(ii,jj) == 1
                    newlast = 1;
                    totalConc = totalConc + 1;
                    
                else
                    
                    newlast = 0;
                end
                
                if (newlast == 1 && last == 0) || (newlast == 0 && last == 1)
                    counter = counter+1;
                end
                last = newlast;
           % end
            
        end
        
    oldmeanFingerWidth = meanFingerWidth;
    meanFingerWidth =2*totalConc/(counter); %* (0.0345);
    meanFingerWidth = max(meanFingerWidth,oldmeanFingerWidth);
    jj = jj + 1;
  end

end

