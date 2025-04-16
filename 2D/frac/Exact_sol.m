function [U] = Exact_sol( X, Y, cx, cy, t, gamma )

sigma = 0.07;
sigma2 = sigma^2;

x0 = 0.5;
y0 = 0.5;

U =(sigma2/(sigma2 + 2*gamma*t)) ...  
     * exp(  - 1/2*((X-cx*t-x0).^2 +  (Y-cy*t-y0).^2)/(sigma2+2*t*gamma));


end

