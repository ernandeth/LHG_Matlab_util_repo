clc
clear all



% Circular wire
% angle = linspace(0,360,20);
% x= 0.020 * cosd(angle);
% y= 0.020 * sind(angle);
% 
% wire(:,1)=x;
% wire(:,2)=y;
% wire(:,3)=0;


% %Figure8coil
angle = linspace(0,360,100);
x1= 0.040* cosd(angle) - 0.040;
y1= 0.040 * sind(angle);


x = [x1(1:end) -x1(1:end)];
y = [y1(1:end) y1(1:end)];


wire(:,1)=x;
wire(:,2)=y;
wire(:,3)=0;


current=1e-3;

% plot3(wire(:,1),wire(:,2),wire(:,3));

[b,bx,by,bz]= biot3d_mod(current,wire);



babs = sqrt(bx.^2 + by.^2 + bz.^2);

% 
% for i=1:64
% %             quiver(bx(:,:,i), by(:,:,i));
%     imagesc(babs(:,:,i)),colorbar;
% %         imagesc(bz(:,:,i)),colorbar;
%     drawnow;
%     pause;
% end


% plot(wire(:,1),wire(:,2)),axis ([-.032 .032 -.032 0.032]),xlabel('x(m)'),ylabel('y(m)');
% figure;imagesc([-0.032:0.004:0.032],[-0.032:0.004:0.032],babs(:,:,7)),axis xy,xlabel('x(m)'),ylabel('y(m)'),colorbar,colormap(gray)
% figure;plot(wire(:,1),wire(:,2)),axis ([-.128 .128 -.128 0.128]),xlabel('x(m)'),ylabel('y(m)');
% figure;imagesc([-0.128:0.004:0.128],[-0.128:0.004:0.128],babs(:,:,32)),axis xy,xlabel('x(m)'),ylabel('y(m)'),colorbar;