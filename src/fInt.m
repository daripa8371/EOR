function out = fInt(T,fmatrix,para,v)

% evaluating source term f at the vertices of the element triangle
f11=f_func(T.x(1),T.y(1),para,fmatrix);
f12=f_func(T.x(2),T.y(2),para,fmatrix);
f13=f_func(T.x(3),T.y(3),para,fmatrix);

out = trmatrix(T,f11,f12,f13,v);
end



