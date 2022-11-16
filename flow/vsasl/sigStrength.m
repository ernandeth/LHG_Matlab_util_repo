function [st_r, angle, vv] = sigStrength(viewpoint, MODE, angle)

switch nargin
    case {0, 1}  % MODE: 0
        viewpoint = [1, 0, 0];
        MODE = int16(0);
        angle = (0:0.5:90)/180 * pi;
        vv = (0:0.05:(4*pi))'+eps;

    case 2  % MODE: mod(MODE, 2)
        MODE = mod(int16(MODE), 2);
        
        switch MODE
            case 0
                angle = (0:0.5:90)/180 * pi;
                vv = (0:0.05:(4*pi))'+eps;
                
            case 1
                viewpoint = [0, 0, 1];
                angle = 90/180 * pi;
                vv_range = (0:0.05:(4*pi))'+eps;
                gRatio = (0:0.05:(4*pi))+eps;
                
                vv = vv_range * gRatio;
                
        end
        
    case 3
        MODE = mod(int16(MODE), 2);
        
        switch MODE
            case 0
                vv = (0:0.05:(4*pi))'+eps;
                
            case 1
                viewpoint = [0, 0, 1];
                vv_range = (0:0.05:(4*pi))'+eps;
                gRatio = (0:0.05:(4*pi))+eps;
                
                vv = vv_range * gRatio;
        end
end

mz = sin(vv)./vv;
mxy = (1-cos(vv))./vv;

st = abs(mz * sin(angle) + mxy*cos(angle));
ref = st(1, :);
ref = ones(size(vv,1), 1) * ref;
st_r = st - ref;

figure(99)

switch MODE
    case 0
        [XX, YY] = meshgrid(angle, vv);
        mesh(XX, YY, -st_r);
        axis([0, pi/2, 0, 4*pi, -1, 1]);
        grid on; colorbar;

        set(gca, 'YTick', 0:pi/2:4*pi); ylabel('velocity');
        set(gca,'YTickLabel',{'0','1/2*pi','pi','3/2*pi','2pi', '5/2*pi', '3pi', '7/2*pi', '4pi'});
        set(gca, 'XTick', 0:pi/8:pi/2); xlabel('tip angle');
        set(gca,'XTickLabel',{'0','1/8*pi','1/4*pi','3/8*pi','1/2*pi'});
        
    case 1
        [XX, YY] = meshgrid(gRatio, vv_range);
        mesh(XX, YY, -st_r);
        axis([0, 4*3.1416, 0, 4*pi, 0, 1]);
        grid on; colorbar;

        set(gca, 'YTick', 0:pi/2:4*pi); ylabel('velocity');
        set(gca,'YTickLabel',{'0','1/2','1','3/2','2', '5/2', '3', '7/2', '4'});
        xlabel('gRatio');
end

view(viewpoint);

end