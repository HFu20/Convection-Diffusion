function  [Mass,u ] = Solve_u_inter_non( miu_x, miu_y,  lam_x, lam_y, nx , ny,dx, dy,  gamma, nt, dt, X, Y, U0, T0,T)
%----matrix ------------- 
[matrix_Delta_x_2,matrix_Delta_x ,matrix_A1x ,matrix_A2x ] = Matrix_interior(nx-1) ;
[matrix_Delta_y_2,matrix_Delta_y ,matrix_A1y ,matrix_A2y ] = Matrix_interior(nx-1);

H_1x = matrix_A1x - gamma*miu_x/4*matrix_Delta_x_2;
H_1y = matrix_A1y - gamma*miu_y/4*matrix_Delta_y_2;
H_2x = matrix_A1x + 1/4*gamma*miu_x*(matrix_Delta_x_2);
H_2y = matrix_A1y + 1/4*gamma*miu_y*(matrix_Delta_y_2);

u_n = U0(2:end-1, 2:end-1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
% kxi_plus = -1e-12 ;

tt = T0 : dt : T ;
Mass_u = zeros( length( tt ),1 ); % 初始质量也要计算
Mass_S = zeros( length( tt ),1 ); 
kxi =0;
M=10;
% N =0 ;
for i =1 : nx-1
    for j = 1 : ny-1
         Mass_u(1) = Mass_u(1)  + u_n(i,j)*dx*dy  ;
        for k = 1 : nt+1
            Mass_S(1) = Mass_S(1) +  fun_Source( X(i),Y(j), gamma, T0+(k-1)*dt  )*dt*dx*dy ;
            %Mass(1) = Mass(1) + u_n(i,j)*dx*dy ; % 只考虑内点质量（边界为0）
        end
    end
end


%%%%%%%%%%%%%%%%%%%

 for k = 1:nt

    %保留A_1y*unPlusHf1,第2步里直接运用
%     u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*matrix_A1y'+ miu_y *gamma/ 2*matrix_A1x*u_n*matrix_Delta_y_2'  );
% 
%     u_nPlusHf2       = (u_bar_nPlusHf1-miu_y*gamma/4* u_n*matrix_Delta_y_2')/ H_1y' ;

    u_bar_nPlusHf1 =  H_1x\( H_2x*u_n*H_2y'  );

    u_nPlusHf2       = u_bar_nPlusHf1/ H_1y' ;
    %-----------limiter------------------------------------------
%        for j=1:ny
%            u_nPlusHf2(:,j) = fun_limiter2sides_final(nx,u_nPlusHf2(:,j),min(min(U0)),max(max(U0)));
%        end
%         for i=1:nx
%           u_nPlusHf2(i,:) = fun_limiter2sides_final(ny,u_nPlusHf2(i,:)',min(min(U0)),max(max(U0)))';
%         end
    %---------------------------------------------------------------
[  u_nPlusHf4 ] = fun_Mstep_convec( u_nPlusHf2, matrix_A2x, matrix_Delta_x, lam_x, lam_y, matrix_Delta_y, matrix_A2y, dt,X, Y, gamma, T0,  k, kxi ,M ) ;

%     u_nPlusHf3    = u_nPlusHf2 - matrix_A2x\matrix_Delta_x*lam_x* fun_f(u_nPlusHf2)...
%         - lam_y*(fun_g(u_nPlusHf2)) *matrix_Delta_y'/matrix_A2y' ...
%         + dt*fun_Source (X,Y, gamma, T0+(k-1)*dt  )  ;  %g
% 
%     %-----------limiter------------------------------------------
% %        for i=1:nx
% %           u_nPlusHf3(i,:) = fun_limiter2sides_final(ny,u_nPlusHf3(i,:)',min(min(U0)),max(max(U0)))';
% %        end
% %         for j = 1:ny
% %           u_nPlusHf3(:,j) = fun_limiter2sides_final(nx,u_nPlusHf3(:, j ),min(min(U0)),max(max(U0)));
% %         end
%     %---------------------------------------------------------------
% 
%     u_nPlusHf4    = u_nPlusHf2 -  matrix_A2x\matrix_Delta_x*lam_x* 0.5*(fun_f(u_nPlusHf2)+fun_f(u_nPlusHf3) ) ...%f
%         -   lam_y/2* (fun_g(u_nPlusHf3)+fun_g(u_nPlusHf2) ) * matrix_Delta_y' /matrix_A2y' ...
%         + dt/2 * ( fun_Source (X,Y, gamma, T0+ (k-1)*dt  ) + fun_Source (X,Y, gamma, T0+(k)*dt  ) ) ;%g;
    %   -----------limiter------------------------------------------
%       for i=1:nx
%         u_nPlusHf4(i,:) = fun_limiter2sides_final(ny,u_nPlusHf4(i,1:end)',min(min(U0)),max(max(U0)))';
%       end
%         for j = 1:ny
%           u_nPlusHf4(:,j) = fun_limiter2sides_final(nx,u_nPlusHf4(:, j ),min(min(U0)),max(max(U0)));
%         end
    % u_nPlusHf4 = u_nPlusHf4';
    %    for i=1:nx
    %     u_nPlusHf4(:,i) = limiter_2sides_final(ny,u_nPlusHf4(:,i),min(min(U0)),max(max(U0)));
    %    end
    % u_nPlusHf4 = u_nPlusHf4';
    %    ---------------------------------------------------------------
%     u_bar_nPlusHf5  =  H_1x\(H_2x*u_nPlusHf4 *matrix_A1y'+miu_y*gamma/ 2*matrix_A1x*u_nPlusHf4 *matrix_Delta_y_2');
% 
%     u_n  =  (u_bar_nPlusHf5-miu_y/4*gamma*u_nPlusHf4*matrix_Delta_y_2')/ H_1y';
    u_bar_nPlusHf5 =  H_1x\( H_2x*u_nPlusHf4*H_2y'  );

   u_nPlus1      = u_bar_nPlusHf5/ H_1y' ;
    %-----------limiter------------------------------------------
%        for j=1:ny
%            u_n(:,j) = fun_limiter2sides_final(nx,u_n(:,j),min(min(U0)),max(max(U0)));
%        end
%         for i=1:nx
%           u_n(i,:) = fun_limiter2sides_final(ny,u_n(i,:)',min(min(U0)),max(max(U0)))';
%         end

    
  %%--------------------------------------------

    %     % Lagrange乘子 preserving bound and mass

    %
 
    for i = 1: nx-1
        for j = 1: ny-1
     Mass_u(k+1) = Mass_u(k+1) + u_nPlus1 (i,j)*dx*dy ;
            Mass_S(k+1) = Mass_S(k+1) + fun_Source ( X(i),Y(j), gamma, T0+(k)*dt  )*dt*dx*dy ;
        end
    end
  
    u_n  = u_nPlus1 ;
end

  Mass =  Mass_u - Mass_S  ;
%-----------limiter------------------------------------------
%    for j=1:ny
%        u_n(:,j) = limiter_2sides_final(nx,u_n(:,j),min(min(U0)),max(max(U0)));
%    end
%     for i=1:nx
%       u_n(i,:) = limiter_2sides_final(ny,u_n(i,:)',min(min(U0)),max(max(U0)))';
%     end
%---------------------------------------------------------------




u = U0 ;
u(2:end-1, 2:end-1) = u_n ;

end

