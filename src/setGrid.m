function out = setGrid(para, U, L,beta)
m = para.box.m;
n = para.box.n;

out = cell(m+1,n+1);
func = @weak;


for j = 1:(m+1)
	for l = 1:(n+1)
        
		
        if j==1&&l~=1&&l~=(n+1)
        t1 = func(L{j,l},para,beta,[1,0,0]);
		t2 = [0,0,0,0];
		t3 = [0,0,0,0];
		t4 = [0,0,0,0];
		t5 = func(L{j,l-1  },para,beta,[0,0,1]);
		t6 = func(U{j,l-1  },para,beta,[0,1,0]);
        end
        if j==(m+1)&&l~=1&&l~=(n+1)
        t1 = [0,0,0,0];
		t2 = func(U{j-1,  l},para,beta,[0,0,1]);
		t3 = func(L{j-1,  l},para,beta,[0,1,0]);
		t4 = func(U{j-1,  l-1  },para,beta,[1,0,0]);
		t5 = [0,0,0,0];
		t6 = [0,0,0,0];  
        end
        if j~=1&&j~=(m+1)&&l==1
        t1 = func(L{j,l},para,beta,[1,0,0]);
		t2 = func(U{j-1,  l},para,beta,[0,0,1]);
		t3 = func(L{j-1,  l},para,beta,[0,1,0]);
		t4 = [0,0,0,0];
		t5 = [0,0,0,0];
		t6 = [0,0,0,0];
        end
        if j~=1&&j~=(m+1)&&l==(n+1)
        t1 = [0,0,0,0];
		t2 = [0,0,0,0];
		t3 = [0,0,0,0];
		t4 = func(U{j-1,  l-1  },para,beta,[1,0,0]);
		t5 = func(L{j,l-1  },para,beta,[0,0,1]);
		t6 = func(U{j,l-1  },para,beta,[0,1,0]);
        end
        if j==1&&l==1
        t1 = func(L{j,l},para,beta,[1,0,0]);
		t2 = [0,0,0,0];
		t3 = [0,0,0,0];
		t4 = [0,0,0,0];
		t5 = [0,0,0,0];
		t6 = [0,0,0,0];		
        end
        if j==1&&l==(n+1)
        t1 = [0,0,0,0];
		t2 = [0,0,0,0];
		t3 = [0,0,0,0];
		t4 = [0,0,0,0];
		t5 = func(L{j,l-1  },para,beta,[0,0,1]);
		t6 = func(U{j,l-1  },para,beta,[0,1,0]);
        end
        if j==(m+1)&&l==1
        t1 = [0,0,0,0];
		t2 = func(U{j-1,  l},para,beta,[0,0,1]);
		t3 = func(L{j-1,  l},para,beta,[0,1,0]);
		t4 = [0,0,0,0];
		t5 = [0,0,0,0];
		t6 = [0,0,0,0];
        end
        if j==(m+1)&&l==(n+1)
        t1 = [0,0,0,0];
		t2 = [0,0,0,0];
		t3 = [0,0,0,0];
		t4 = func(U{j-1,  l-1  },para,beta,[1,0,0]);
		t5 = [0,0,0,0];
		t6 = [0,0,0,0];
        end
        if j~=1&&j~=(m+1)&&l~=1&&l~=(n+1)
		t1 = func(L{j,l},para,beta,[1,0,0]);
		t2 = func(U{j-1,  l},para,beta,[0,0,1]);
		t3 = func(L{j-1,  l},para,beta,[0,1,0]);
		t4 = func(U{j-1,  l-1  },para,beta,[1,0,0]);
		t5 = func(L{j,l-1  },para,beta,[0,0,1]);
		t6 = func(U{j,l-1  },para,beta,[0,1,0]);
        end
%         %%%% debugging and documenting code
%         t1
%         t2
%         t3
%         t4
%         t5
%         t6
%         %%%%%%
		grid.c = t1(1) + t2(3) + t3(2) + t4(1) + t5(3) + t6(2);
		grid.w = t3(1) + t4(2);
		grid.s = t4(3) + t5(1);
		grid.n = t1(3) + t2(1);
		grid.e = t1(2) + t6(1);
		grid.nw = t2(2) + t3(3);
		grid.se = t5(2) + t6(3);
		grid.const = t1(4) +t2(4) +t3(4) +t4(4) +t5(4) +t6(4);
%         disp(grid)
		out{j,l} = grid;
%         pause
	end
end