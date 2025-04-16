function [S] = fun_Source( x, y, gamma , t,sigma2 )
% sigma2 = 0.001/2;
x0 = -0.35 ;
y0 = 0  ;
%  S = -2..*pi..*( (x-cos(pi/2..*y)..*sin(pi/2..*x))..*( x0..*sin(pi..*t) + y0..*cos(pi..*t) -y ) + ...
%                      (y-cos(pi/2..*x)..*sin(pi/2..*y))..*(-x0..*cos(pi..*t) + y0..*sin(pi..*t)+x))..*u ;
[Y,X] = meshgrid(x,y) ;
S = -(sigma2.*pi.*exp(-((y0 - Y.*cos(pi.*t) + X.*sin(pi.*t)).^2 + (X.*cos(pi.*t) - x0 + Y.*sin(pi.*t)).^2) ...
    /(2.*sigma2 + 4.*gamma.*t)).*(X.*y0.*cos(pi.*t) - x0.*Y.*cos(pi.*t) + X.*x0.*sin(pi.*t) + Y.*y0.*sin(pi.*t)  ...
    - X.*cos((pi.*X)/2).*sin((pi.*Y)/2) + Y.*cos((pi.*Y)/2).*sin((pi.*X)/2) + x0.*cos(pi.*t).*cos((pi.*X)/2).*sin((pi.*Y)/2)  ...
    - y0.*cos(pi.*t).*cos((pi.*Y)/2).*sin((pi.*X)/2) - x0.*cos((pi.*Y)/2).*sin(pi.*t).*sin((pi.*X)/2)  ...
    - y0.*cos((pi.*X)/2).*sin(pi.*t).*sin((pi.*Y)/2)))/(sigma2 + 2.*gamma.*t).^2 ;
 

end

 

