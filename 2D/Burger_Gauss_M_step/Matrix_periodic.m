function [matrix_Delta_x_2,matrix_Delta_x ,matrix_A1 ,matrix_A2 ] = Matrix_periodic(nx)

%周期型
matrix_Delta_x_2 =  sparse(1:nx, 1:nx, -2, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1, nx, nx);
matrix_Delta_x_2(1,end)=1;
matrix_Delta_x_2(end,1)=1;

matrix_Delta_x = sparse(1:nx-1, 2:nx, 1/2, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, -1/2, nx, nx);
matrix_Delta_x(1,end)=-1/2;
matrix_Delta_x(end,1)=1/2;
              
%matrix_I = sparse(1:nx, 1:nx, 1, nx, nx);
matrix_A1 =sparse(1:nx, 1:nx, 5/6, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1/12, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1/12, nx, nx);
matrix_A1(1,end)=1/12;
matrix_A1(end,1)=1/12;

matrix_A2 =sparse(1:nx, 1:nx, 2/3, nx, nx) ... 
                  + sparse(1:nx-1, 2:nx, 1/6, nx, nx) ... 
                  + sparse(2:nx, 1:nx-1, 1/6, nx, nx);
matrix_A2(1,end)=1/6;
matrix_A2(end,1)=1/6;



end

