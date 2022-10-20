function out = weak(T,para,beta,v)
    % evaluating beta at the vertices of the element triangle
    beta1=beta_func(T.x(1),T.y(1),para,beta);

    beta2=beta_func(T.x(2),T.y(2),para,beta);

    beta3=beta_func(T.x(3),T.y(3),para,beta);
    
    % computing average of the beta values at the vertices
    beta10=(beta1+beta2+beta3)/3;
    out = weak1(T,beta10,v);
   
end