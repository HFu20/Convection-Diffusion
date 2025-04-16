function   [Mass,u, kxi_plus,yinta, Gplus, Count] = Solve_u_inter_DYak_S( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X, Y, U0, T0,T,K)

%----matrix -------------
[ matrix_Delta_x_2,matrix_Delta_x ,matrix_A1x ,matrix_A2x ] = Matrix_interior( nx-1 );
[ matrix_Delta_y_2,matrix_Delta_y ,matrix_A1y ,matrix_A2y ] = Matrix_interior( ny-1 );

H_1x = matrix_A1x - gamma*miu_x/4*matrix_Delta_x_2;
H_1y = matrix_A1y - gamma*miu_y/4*matrix_Delta_y_2;
H_2x = matrix_A1x + 1/4*gamma*miu_x*(matrix_Delta_x_2);
H_2y = matrix_A1y + 1/4*gamma*miu_y*(matrix_Delta_y_2);

u_n = U0( 2:end-1, 2:end-1 ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
% kxi_plus = -1e-12 ;
kxi = 0 ;
tt = T0 : dt : T ;
Mass = zeros( length( tt ),1 ); % 初始质量也要计算
Count = zeros( nt ,1);
kxi_plus = zeros( nt ,1 );
Gplus = zeros( nt ,1 );
Gminus = zeros( nt ,1 );
yinta = zeros( nt ,1 );

% N =0 ;
% for i =1 : nx -1
%     for j = 1 : ny-1
%         Mass(1) = Mass(1) + u_n(i,j)*dx*dy ; % 只考虑内点质量（边界为0）
%     end
% end
%

%%%%%%%%%%%%%%%%%%%

tic
 for k = 1:nt

    u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*H_2y'  );

    u_nPlusHf2       =   u_bar_nPlusHf1/ H_1y' ;
  
    [  u_nPlusHf4 ]  =   fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k, kxi ,K ) ;

    u_bar_nPlusHf5 =   H_1x\( H_2x*u_nPlusHf4*H_2y'  );

   u_nPlus1              =    u_bar_nPlusHf5/ H_1y' ;

    %     % Lagrange乘子 preserving bound and mass

    %
%     [ kxi_plus( k ), Gplus( k ),   Count(k ) ] = fun_kxi_plus_utidle( dt,  u_nPlus1, u_n, min(min(U0)), max(max(U0)),  kxi, dx,dy );
%     yinta(k) = dt*( kxi_plus(k) - kxi ) ;
%     z = u_nPlus1+ yinta(k) ;
%     kxi = kxi_plus( k ) ;
%     for i = 1: nx-1
%         for j = 1: ny-1
% 
% %             if z( i,j ) <= min( min(U0) )
% %                 u_nPlus1 ( i,j ) = min(min(U0));
% %             elseif z( i,j ) >= max(max(U0))
% %                 u_nPlus1 (i,j) = max(max(U0));
% %             else
% %                 u_nPlus1 (i,j) = z(i, j ) ;
% %             end
% 
%             Mass(k+1) = Mass(k+1) + u_nPlus1 (i,j)*dx*dy ;
%         end
%     end

    
    u_n  = u_nPlus1 ;
end


%-----------limiter------------------------------------------
%    for j=1:ny
%        u_n(:,j) = limiter_2sides_final(nx,u_n(:,j),min(min(U0)),max(max(U0)));
%    end
%     for i=1:nx
%       u_n(i,:) = limiter_2sides_final(ny,u_n(i,:)',min(min(U0)),max(max(U0)))';
%     end
%---------------------------------------------------------------



u = U0 ;
u(2:end-1,2:end-1) = u_n ;



end


