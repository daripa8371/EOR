function B = setB(para,grid,rh)
m = para.box.m;
n = para.box.n;

B = zeros((m+1)*(n+1),1);

for j = 1 : (m+1)
	for l = 1 : (n+1)
		a = grid{j,l};
		id = j+(l-1)*(m+1);
		B(id) = a.const;
% 		if j == 1 && l ~= 1 && l ~= (n-1)	
% 			B(id) = B(id) + a.w * value(para,0,l,phi) + a.nw * value(para,0,l+1,phi);
% 		end
% 		if j == 1 && l == 1
% 			B(id) = B(id) + a.w*value(para,0,1,phi) + a.nw*value(para,0,2,phi) + a.s*value(para,1,0,phi) + a.se*value(para,2,0,phi);
% 		end
% 		if j == 1 && l == (n-1)
% 			B(id) = B(id) + a.w * value(para,0,n-1,phi) + a.nw * value(para,0,n,phi) + a.n * value(para,1,n,phi);
% 		end
% 		if j == (m-1) && l ~= 1 && l ~= (n-1)	
% 			B(id) = B(id) + a.e * value(para,m,l,phi) + a.se * value(para,m,l-1,phi);
% 		end
% 		if j == (m-1) && l == 1
% 			B(id) = B(id) + a.e * value(para,m,1,phi) + a.se * value(para,m,0,phi) + a.s * value(para,m-1,0,phi);
% 		end
% 		if j == (m-1) && l == (n-1)
% 			B(id) = B(id) + a.e*value(para,m,n-1,phi) + a.se*value(para,m,n-2,phi) + a.n*value(para,m-1,n,phi) + a.nw*value(para,m-2,n,phi);
% 		end
% 		if l == 1 && j ~= 1 && j~= (m-1)
% 			B(id) = B(id) + a.s * value(para,j,0,phi) + a.se * value(para,j+1,0,phi);
% 		end
% 		if l == (n-1) && j ~= 1 && j~= (m-1)
% 			B(id) = B(id) + a.n * value(para,j,n,phi) + a.nw * value(para,j-1,n,phi);
% 		end
	end
end


B = rh - B;
