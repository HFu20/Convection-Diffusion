function [U] = Exact_sol(x,t,gamma,a,c)
%-----------------------------------------------------------------
%exact solution:   b=gamma;
%----------------------------------------------------------------------------
U = a-c*tanh(c/(2*gamma)*(x-a*t)) ; % Initial /Finitial时刻精确解
end










