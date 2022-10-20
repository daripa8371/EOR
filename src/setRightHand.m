function out = setRightHand(para,fmatrix,U, L)

m = para.box.m;
n = para.box.n;
	
rh = zeros((m+1)*(n+1),1);
func = @fInt;
% sum=0;
for j = 1:(m+1)
	for l = 1:(n+1)
		id = j + (l-1)*(m+1);        
        if j==1&&l~=1&&l~=(n+1)
        t1 = func(L{j,l},fmatrix,para,[1,0,0]);
		t2 = 0;
		t3 = 0;
		t4 = 0;
		t5 = func(L{j,l-1  },fmatrix,para,[0,0,1]);
		t6 = func(U{j,l-1  },fmatrix,para,[0,1,0]);
        end
        if j==(m+1)&&l~=1&&l~=(n+1)
        t1 = 0;
		t2 = func(U{j-1,  l},fmatrix,para,[0,0,1]);
		t3 = func(L{j-1,  l},fmatrix,para,[0,1,0]);
		t4 = func(U{j-1,  l-1  },fmatrix,para,[1,0,0]);
		t5 = 0;
		t6 = 0;  
        end
        if j~=1&&j~=(m+1)&&l==1
        t1 = func(L{j,l},fmatrix,para,[1,0,0]);
		t2 = func(U{j-1,  l},fmatrix,para,[0,0,1]);
		t3 = func(L{j-1,  l},fmatrix,para,[0,1,0]);
		t4 = 0;
		t5 = 0;
		t6 = 0;
        end
        if j~=1&&j~=(m+1)&&l==(n+1)
        t1 = 0;
		t2 = 0;
		t3 = 0;
		t4 = func(U{j-1,  l-1  },fmatrix,para,[1,0,0]);
		t5 = func(L{j,l-1  },fmatrix,para,[0,0,1]);
		t6 = func(U{j,l-1  },fmatrix,para,[0,1,0]);
        end
        if j==1&&l==1
        t1 = func(L{j,l},fmatrix,para,[1,0,0]);
		t2 = 0;
		t3 = 0;
		t4 = 0;
		t5 = 0;
		t6 = 0;
        end
        if j==1&&l==(n+1)
        t1 = 0;
		t2 = 0;
		t3 = 0;
		t4 = 0;
		t5 = func(L{j,l-1  },fmatrix,para,[0,0,1]);
		t6 = func(U{j,l-1  },fmatrix,para,[0,1,0]);
        end
        if j==(m+1)&&l==1
        t1 = 0;
		t2 = func(U{j-1,  l},fmatrix,para,[0,0,1]);
		t3 = func(L{j-1,  l},fmatrix,para,[0,1,0]);
		t4 = 0;
		t5 = 0;
		t6 = 0;
        end
        if j==(m+1)&&l==(n+1)
        t1 = 0;
		t2 = 0;
		t3 = 0;
		t4 = func(U{j-1,  l-1  },fmatrix,para,[1,0,0]);
		t5 = 0;
		t6 = 0;
        end
        if j~=1&&j~=(m+1)&&l~=1&&l~=(n+1)
		t1 = func(L{j,l},fmatrix,para,[1,0,0]);
		t2 = func(U{j-1,  l},fmatrix,para,[0,0,1]);
		t3 = func(L{j-1,  l},fmatrix,para,[0,1,0]);
		t4 = func(U{j-1,  l-1  },fmatrix,para,[1,0,0]);
		t5 = func(L{j,l-1  },fmatrix,para,[0,0,1]);
		t6 = func(U{j,l-1  },fmatrix,para,[0,1,0]);
        end
        
		rh(id) = t1 + t2 + t3 + t4 + t5 + t6;
	end
end
% sum
out = rh;


