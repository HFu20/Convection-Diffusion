function [matrix_LHS_step1, matrix_RHS_step1, matrix_A2, matrix_Delta_x] = Matrix_nonp(nx ,dt, dx, gamma)
 
matrix_Delta_x_2 =  sparse(1:nx, 1:nx, -2, nx, nx) ...
    + sparse(1:nx-1, 2:nx, 1, nx, nx) ...
    + sparse(2:nx, 1:nx-1, 1, nx, nx);
matrix_Delta_x_2 = matrix_Delta_x_2 /dx^2 ;

matrix_Delta_x = sparse(1:nx-1, 2:nx, 1/2, nx, nx) ...
                              + sparse(2:nx, 1:nx-1, -1/2, nx, nx);
matrix_Delta_x = matrix_Delta_x/dx ;

%matrix_I = sparse(1:nx+1, 1:nx+1, 1, nx+1, nx+1);
matrix_A1 =sparse(1:nx, 1:nx, 5/6, nx, nx) + sparse(1:nx-1, 2:nx, 1/12, nx, nx) ...
                     + sparse(2:nx, 1:nx-1, 1/12, nx, nx);


matrix_A2 = sparse(1:nx, 1:nx, 2/3, nx, nx)  + sparse(1:nx-1, 2:nx, 1/6, nx, nx) ...
                      + sparse(2:nx, 1:nx-1, 1/6, nx, nx);
 


matrix_LHS_step1 = matrix_A1 - dt/2 * gamma/2  * matrix_Delta_x_2;
matrix_RHS_step1 = matrix_A1 + dt/2 * gamma/2  * matrix_Delta_x_2;

end

