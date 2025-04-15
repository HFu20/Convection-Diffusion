function [S] = fun_Source_wan( x,t )
%% when  gamma = b = 1/64
  %S = zeros( length(x),1);
%% when  gamma = b = 1/10000
S =  100*tanh(750*t - 500*x).*(tanh(750*t - 500*x).^2 - 1) ;
end

