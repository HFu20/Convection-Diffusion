function [ f ] = fun_f( u ,x,y)
% linear
%  f =  u ;
 % nonlinear 
%    f = 0.5 .* u.* u ; 
%     f = sin(pi*y).*u ;   %效果与sin(pi*ones(size(y))'*y).*u一致，但需注意【Y1,X1】
%  f = sin(pi*ones(size(y))'*y).*u ; 

f = -pi*cos((pi*x')/2)*sin((pi*y)/2).*u ;

end

