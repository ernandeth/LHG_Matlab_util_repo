%for mysub = {'RH' 'AB' 'DC' 'PD' 'PS' 'KH' 'MS' 'NS' 'RM' 'HP' 'LB'}
%for mysub = {'RH' 'TD'}

clear mysub
mysub{1} = 'biman5_sb';
%img = readim2([mysub{1} '3dvol']);
img = readim2(['t1spgr']);

fig=figure;
set(fig,'DoubleBuffer','on');
set(gcf,'Color','k')
colormap(gray)
set(gca,'NextPlot','replace','Visible','off')
axis off

%try
	clear mov
	mov = avifile([mysub{1} '.avi'],'Quality',50,'Compression','Indeo3','Fps',3)
     
	for i=1:size(img,3)
		imagesc(img(:,:,i)); axis image;
		drawnow
       	% set(h,'EraseMode','xor');
       	F = getframe(gca);
       	mov = addframe(mov,F);
	end
       
	mov = close(mov);
	
    %catch
	
	%mov = close(mov);

    %end




%end
