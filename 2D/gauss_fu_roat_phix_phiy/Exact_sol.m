function [u] = Exact_sol( x, y,  t, gamma,sigma2  )


% sigma2 = 0.005/2;

% x0 = -0.4 ;
x0 = -0.35 ;
y0 = 0  ;

xbar = x'*ones(size(y))*cos(pi*t) + ones(size(x'))*y*sin(pi*t) ;
ybar = -x'*ones(size(y))*sin(pi*t) + ones(size(x'))*y*cos(pi*t) ;
u = 2*sigma2/(  2*sigma2 + 4*gamma*t )*exp(  -(  (xbar-x0).^2  +  (ybar-y0).^2 )/( 2*sigma2 + 4*gamma*t) );
 
end

