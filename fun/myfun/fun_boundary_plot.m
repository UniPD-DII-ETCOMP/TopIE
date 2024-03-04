function [none] = fun_boundary_plot(F1_m,F1_e,Matrix_P0_m,Matrix_P0_e,alpha,...
                                    id_fill,id_empty,id_figure)
gray=[0 0 0]; %[0.43 0.43 0.43];
green=[0 0.3 0];
figure(id_figure)
hold on
grid on
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor',gray,'FaceAlpha',alpha) 
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor',green,'FaceAlpha',0.25) 
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor','r','FaceAlpha',alpha) 
patch('Faces',F1_m(1:4,id_fill).','Vertices',Matrix_P0_m,'Facecolor',gray,'FaceAlpha',alpha,'edgecolor','none') 
patch('Faces',F1_m(1:4,id_empty).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',0.25,'edgecolor','none') 
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'Facecolor','r','FaceAlpha',alpha) 
legend('filled','empty','winding')
axis equal
title('Mesh')
xlabel('r (m)')
ylabel('\phi')
zlabel('z (m)')
set(gca,'FontName','Century','FontSize',16);
view(0,0)
Matrix_P0=[Matrix_P0_m;Matrix_P0_e];
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) ...
    min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 ...
    min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
drawnow
end