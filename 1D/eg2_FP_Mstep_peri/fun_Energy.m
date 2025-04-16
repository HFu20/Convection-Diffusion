function [E] = fun_Energy( x,u, dx )
 
E = ( 0.5*x.*x.*u  + u.*log(u) + (1-u).*log( 1-u ) ) ;
E = sum( E* dx ) ;

end

