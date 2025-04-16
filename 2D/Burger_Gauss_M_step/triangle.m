 function triangle(tr1,tr2,ratio)
% tr1：水平边和斜边交点坐标
% tr2：水平边另一个端点坐标
% ratio：斜边斜率
% text1：水平边文字
% text2：垂直边文字
% pos adjust：文字坐标的微调   
    x0 = tr1(1); y0 = tr1(2); 
    x1 = tr2;  y1 = y0*(x1/x0)^ratio;
%不希望这样绘制的线条在添加图例时被探测到
    loglog([x0,x1],[y0,y0],'-k','HandleVisibility','off','Linewidth',1);
    loglog([x1,x1],[y0,y1],'-k','HandleVisibility','off','Linewidth',1);
    loglog([x1,x0],[y1,y0],'-k','HandleVisibility','on','Linewidth',1);
%     text(sqrt(x0*x1),y0*pos_adjust);
%     text(x1*pos_adjust,sqrt(y0*y1));
end