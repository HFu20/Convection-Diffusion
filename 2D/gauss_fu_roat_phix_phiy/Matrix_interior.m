function [matrix_Delta_x_2,matrix_Delta_x ,matrix_A12 ,matrix_A6 ] = Matrix_interior(nx)

%周期型
matrix_Delta_x_2 =  sparse(1:nx, 1:nx, -2, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1, nx, nx);
% matrix_Delta_x_2=matrix_Delta_x_2/dx^2;

matrix_Delta_x = sparse(1:nx-1, 2:nx, 1/2, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, -1/2, nx, nx);
% matrix_Delta_x = matrix_Delta_x/(dx) ;
             
%matrix_I = sparse(1:nx, 1:nx, 1, nx, nx);
matrix_A12 =sparse(1:nx, 1:nx, 5/6, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1/12, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1/12, nx, nx);


matrix_A6 =sparse(1:nx, 1:nx, 2/3, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1/6, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1/6, nx, nx);



end

