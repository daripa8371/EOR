function A = setA(para,grid)
m = para.box.m;
n = para.box.n;

A = zeros((m+1)*(n+1)*7,3);
list = 0;
for j = 1 : (m+1)
	for l = 1 : (n+1);
		a = grid{j,l};
		id = j+(l-1)*(m+1);
		list = list+1;
		A(list,:) = [id,id,a.c];
		if j ~= 1
			list = list+1;
			A(list,:) = [id,id-1,a.w];
		end
		if j ~= 1 && l ~= (n+1)
			list = list+1;
			A(list,:) = [id,id+m,a.nw];
		end
		if l ~= (n+1)
			list = list+1;
			A(list,:) = [id,id+m+1,a.n];
		end
		if j ~= (m+1)
			list = list+1;
			A(list,:) = [id,id+1,a.e];
		end
		if l ~= 1
			list = list+1;
			A(list,:) = [id,id-m-1,a.s];
		end
		if j ~= (m+1) && l ~= 1
			list = list+1;
			A(list,:) = [id,id-m,a.se];
		end
	end
end

A = A(1:list,:);