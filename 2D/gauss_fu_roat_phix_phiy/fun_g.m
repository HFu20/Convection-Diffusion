function [ g ] = fun_g( u, x,y )
% linear 
%  g = u ;
% nonliear 
%   g =  0.5.* u .*u ;
%   g = sin(pi*x).*u ;
%     g = sin(pi*x).*u ;

g = pi*sin((pi*x')/2)*cos((pi*y)/2).*u;
end

