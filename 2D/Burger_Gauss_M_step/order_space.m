
%求空间阶数
clear
%------------------
   N_x = [ 60 120  240 360 ];
%  N_x = [   50 100 150 200];
% gamma = 1/64;
% N_x = [ 200, 250,300,350];
%  N_x = [ 150,200, 250,300];
% N_x = [ 90,180,270,360];
 gamma = 0.005;
df = 1 ;
cx = 1 ;  cy = 1 ;

xmin = -1 ; xmax = 2 ;
ymin = -1 ; ymax = 2;

T0 = 0 ; T_old   = 0.6 ;
%%
nstep_x = length( N_x ) ;
Error_Table = zeros( nstep_x ,5 );

for i = 1: nstep_x
    
    nx = N_x(i);
    ny = N_x(i);
  
dx = ( xmax - xmin )/ nx ;
dy = ( ymax - ymin ) / ny ;

X = xmin : dx : xmin + nx*dx ;
Y = ymin : dy : ymin + ny*dy ;
x = X( 1: end-1 ) ; 
y = Y(1: end-1 ) ;
% [ X1, Y1 ] = meshgrid( Y , X ) ; % (nx+1) * (nx+1)
[ X2, Y2 ] = meshgrid( y, x ) ; %  

%  dt = dx/( 6*df ) ;
%   dt =dx*dx;
nt = nx^2 ;
% 
dt =   ( T_old - T0 )/ nt   ;
 T = T0 + nt*dt ;
if T ~= T_old
    fprintf(' Wrong! T_rel = %f \n ', T) 
end
tt = T0 : dt : T ;
%-- Initial ------------
U0 = Exact_sol( X2, Y2, cx, cy, T0, gamma ) ;
U  =  Exact_sol( X2, Y2, cx, cy, T, gamma ) ;
% u_n = U0  ;

%% common part
miu_x=dt/(dx*dx);
miu_y=dt/(dy*dy);
lam_x=dt/dx;
lam_y=dt/dy;
M=1;
% [u] = fun_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny, gamma, nt, dt, X2, Y2, U0, T0);
tic
% [Mass_u,u, kxi_plus,yinta, Gplus, Mass_err_u] = mass_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X2, Y2,x,y, U0, T0,M) ;
 [u, Mass_err_u ] = u_Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,   gamma, nt, dt, X2, Y2,  U0, T0, M,T,dx,dy) ;
v1(i) = toc ;
%-----------------------------------
% figure   
% mesh( X1, Y1, U ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Exact solution")
% 
%  figure 
% mesh( X1, Y1, u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Numerical solution")
% 
%  figure 
% mesh( X1, Y1, U-u ) ;
% xlabel("x-axis"); ylabel("y-axis")
% title("Error")
% 
 figure 
plot( tt(2:end), Mass_err_u ) ;
xlabel("x-axis"); ylabel("y-axis")
title("Error")
%-----------------------
Emax = max(max(abs(U-u))) ;
E2 = sqrt(sum(dy*(sum(dx*((U-u).^2))))) ;


    % put n in the 1st column of table
    Error_Table(i,  1) = nx;
    % compute the Lmax norm of error  and put it in the 4th column
    Error_Table(i, 2)  =  Emax;
    % compute the L2 norm of error and put it in the 2nd column
    Error_Table(i, 4)  =  E2 ;
 

    if i ~=1
        % compute the ratio for Emax and put it in the 5th column of table
        Error_Table(i, 3) = -log(Error_Table(i, 2)/Error_Table(i-1, 2))...
            /log(Error_Table(i, 1)/Error_Table(i-1, 1));

        % compute the ratio for E2, and put it in the 5th column of table
        Error_Table(i, 5) = -log(Error_Table(i, 4)/Error_Table(i-1, 4))...
            /log(Error_Table(i,1)/Error_Table(i-1, 1));

    end


end

format short e
% convert array to table for visualization
TableResults_Err = array2table(Error_Table, 'VariableNames', ...
     {'nx', 'E_max', 'ratio E_max','E_2','ratio E_2'})
%--------------------------------------------------
 
fprintf('All time = %f \n',v1 )
writetable( TableResults_Err,'order_log.xls' ) ;
% figure(3)
% plot(log10(N_x),log10(Error_Table(:, 2)),'b-.^','linewidth',1);
% hold on
% plot(log10(N_x),log10(Error_Table(:, 4)),'ro-','linewidth',1);
% xlabel('log_{10}(N)');
% ylabel('log_{10}(error)');
% % xx=[6 6 6.5 6];
% % yy=[-14 -16 -16 -14];
% % line(xx,yy);
% lgd=legend('L_{\inf}','L_{2}');
% lgd.FontSize = 12;
% set(gca,'FontSize',12);
% % set(gca,'linewidth',0.7);




figure; 
%set(gcf,'unit','centimeters','position',[10,10,20,20]);
%  xx =  logspace(N_x(1), N_x(end), nstep_x);
% loglog(xx, Error_Table(:, 2),'b-o','LineWidth',1.2);
% hold on;
% loglog(xx, Error_Table(:, 4),'r-.^','LineWidth',1.2);
figure
xx = N_x ;
semilogy( xx,  Error_Table(:,2), 'b-o',LineWidth=1.2) ;
hold on 
semilogy( xx, Error_Table(:,4),'r-.^',LineWidth=1.2 ) ;
xlabel('N')
ylabel('error')
tr1 = [150 , 2e-6 ]; tr2 =95 ; ratio = -4;
triangle(tr1,tr2,ratio);
 legend('L_{inf} error','L_{2} error','{Ratio =-4}');
% figure
% xx=[55 80 55 55];
% yy=[5e-6 5e-6 5e-4 5e-6];
% line(xx,yy);

% fid = fopen('Order_loglog.txt', 'w');
% 
% % 将矩阵写入文件
 
% fid = fopen('order_log.txt','wt');
% fprintf(fid,'%r\n',TableResults_Err);      
% fclose(fid);
 
