function out = eval_bs(Q,C,miuw,miuo,flag)
if flag == 0 % 1st formulation of f(s,c)
    out = ((3*miuo*miuw*Q.^2).*((1-Q).^2))./(((miuo*Q.^3) + (miuw*(1-Q).^3)).^2);
else  % 2nd formulation of f(s,c)
    out = ((3*(0.5+C).*Q.^2).*((1-Q).^2))./(((Q.^3) + ((0.5+C).*(1-Q).^3)).^2);
end

% i = 1:size(Q,2); j = 1:size(Q,1);
% out = (3*miuo*miuw*Q(j,i)^2)*((1-Q(j,i))^2)/(((miuo*Q(j,i)^3) + (miuw*(1-Q(j,i))^3))^2);