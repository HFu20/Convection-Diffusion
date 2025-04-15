function [  matrix_Hl , matrix_Hr ,matrix_Delta_x_2 , matrix_Delta_x, matrix_A1,matrix_A2] = Matrix_nonp(n, gamma, dx, dt)

matrix_Delta_x_2 =  sparse(1:n, 1:n, -2, n, n) ...
    + sparse(1:n-1, 2:n, 1, n, n) ...
    + sparse(2:n, 1:n-1, 1, n, n);
matrix_Delta_x_2  = matrix_Delta_x_2 /dx^2 ;


matrix_Delta_x = sparse(1:n-1, 2:n, 1/2, n, n) ...
                              + sparse(2:n, 1:n-1, -1/2, n, n);
matrix_Delta_x = matrix_Delta_x /dx ;

%matrix_I = sparse(1:nx+1, 1:nx+1, 1, nx+1, nx+1);
matrix_A1 =sparse(1:n, 1:n, 5/6, n, n) + sparse(1:n-1, 2:n, 1/12, n, n) ...
                     + sparse(2:n, 1:n-1, 1/12, n, n);


matrix_A2 = sparse(1:n, 1:n, 2/3, n, n)  + sparse(1:n-1, 2:n, 1/6, n, n) ...
                      + sparse(2:n, 1:n-1, 1/6, n, n);

 matrix_Hl = matrix_A1 -  gamma*dt /4 * matrix_Delta_x_2;
 matrix_Hr = matrix_A1 + gamma*dt /4  * matrix_Delta_x_2;


end

