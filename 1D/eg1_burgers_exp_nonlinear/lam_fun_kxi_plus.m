function [ kxi_plus,count ] = lam_fun_kxi_plus( dt,  u_nPlus1,u_n,U0, kxi, dx,  lam_n  )
%额外加入了lambda
kxi_minus = 0; % 赋初值
kxi_k = dt ; %赋初值
w = dx ;  %权重
% eps =  1.0e-20; %容许误差
count = 0; %迭代次数计数器


a = min(U0) ;
b = max(U0) ;
g = @(x) (b-x).*(x-a) ;
%----------------------------------
% 先不看时间层的变化，都是在第 n+1 层迭代，
% 以通过第 n 层的 kxi 计算出新的kxi_plus (n+1层)
%-------------------------------------------
yinta_minus = 1/1*dt*(   kxi_minus - kxi -  lam_n.*g(u_n)   ) ; % alp_k = 1;
v_minus = u_n+  yinta_minus ; %  原文多了dt
Gminus_1 = 0 ; %表示( k-1 )层的四个求和项（只需计算一遍）
Gminus_2 = 0 ;
Gminus_3 = 0 ;
Gminus_4 = 0;

for i = 1: length( u_nPlus1 )

    if  v_minus(i) <= a
        Gminus_3 = Gminus_3+ a* w ;
    elseif  v_minus(i) >= b
        Gminus_2 = Gminus_2 +b* w;
    else
        Gminus_1 = Gminus_1 + v_minus(i)* w ;
    end

    Gminus_4 = Gminus_4 + u_n(i)* w ;

end
Gminus = Gminus_1 + Gminus_2 + Gminus_3 - Gminus_4;

while abs( kxi_k - kxi_minus ) > 1e-16  %eps =  1.0e-6;

    Gk_1 = 0 ;   %表示第k层的四个求和项，需复制后迭代计算
    Gk_2 = 0 ;
    Gk_3 = 0 ;
    Gk_4 = 0;

    count = count + 1;
    if count > 20
        disp('error:不收敛');
        return;
    end

    % 需要用到 Fn( kxi_k ) = Gn( yinta_k )和 Fn( kxi_minus ) = Gn( yinta_minus )，
    % 即 Gk 和 Gminus
    yinta_k = 1/1*dt*(   kxi_k - kxi -lam_n.*g(u_n) ) ; %对于时间二阶格式，B_{k-1} = 1 ;
    v_k = u_nPlus1 + yinta_k ;
    % yinta_k 为判断所在哪个区域的自变量，v_k为判断的整体式子及求和的值
    % 给出 Gk 和 Gminus 的表达式
    for i = 1: length( u_nPlus1 )
        if  v_k(i) <= a
            Gk_3 = Gk_3+ a *w ; % w =1;
        elseif  v_k(i) >= b
            Gk_2 = Gk_2 +b* w ;
        else
            Gk_1 = Gk_1 + v_k(i)* w ;
        end
        Gk_4 = Gk_4 + u_n( i)* w ;  % （u^n,1）的内积
    end
    Gk = Gk_1 + Gk_2 + Gk_3 - Gk_4;

    %     迭代求解 kxi
    if Gk - Gminus  ~= 0
         kxi_plus = kxi_k - Gk*( kxi_k - kxi_minus ) / ( Gk - Gminus )  ;  % 分母=0
%          kxi_plus = kxi_k + eps*sign( kxi_k - kxi_minus )/5;
         %     disp( Gk - Gminus )
    else
    %     kxi_plus = kxi_k + eps*sign( kxi_k - kxi_minus )/5; %效果也可以
        kxi_plus = kxi_k  ;
        fprintf( 'Gk = Gminus at %d *dt \n ',count )
        %continue
        break
    end

    kxi_minus = kxi_k ;
    kxi_k = kxi_plus;
    Gminus =  Gk;




end

% X = [name,' will be ',num2str(age),' this year.'];
% disp(['count is %d',count])
% fprintf( 'count is %d.\n ',count )

end


