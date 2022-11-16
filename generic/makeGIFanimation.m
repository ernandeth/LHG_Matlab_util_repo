if loopidx == 1
    set(gcf,'WindowState','fullscreen')
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if itr == 1
        imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',1/framerate);
    else
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',1/framerate);
    end
end
itr = itr+1;