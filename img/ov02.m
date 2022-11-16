function [fig1, fig2, fig3] =  ov02(h,d,x,y,z,roi, cmin, cmax)
%function [fig1, fig2, fig3] =  ov02(h,d,x,y,z,roi)
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

if isempty(h)
    fig1=subplot (223);
    imagesc(squeeze(d(:,:,z)))  ,axis xy %,axis tight
    caxis([cmin cmax])
    hold on; plot(y,x,'go');hold off;
    set(gca,'Position',[ 0.05    0.01    0.45    0.45]);

    fig2=subplot (222);
    imagesc(squeeze(d(:,y,:))'),axis xy %, axis tight
    caxis([cmin cmax])
    hold on; plot(x,z,'go');hold off;
    set(gca,'Position',[ 0.5    0.5    0.45    0.45]);

    fig3=subplot (221);
    imagesc(squeeze(d(x,:,:))'), axis xy%, axis tight
    caxis([cmin cmax])
    hold on; plot(y,z,'go');hold off;
    set(gca,'Position',[ 0.05    0.5    0.45    0.45]);
    
    if roi>0
        value=mean(mean(d(x-roi:x+roi, y-roi:y+roi,z-roi:z+roi)));
    else
        value=d(x,y,z);
    end
    fprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, value);
else
    stretch = abs(h.zsize/h.xsize);
    %colordef black

    fig1=subplot (223);
    imagesc(squeeze(d(:,:,z))), axis ([0.5 h.ydim 0.5 h.xdim]) , axis xy %xy %,axis tight
    caxis([cmin cmax])
    hold on; plot(y,x,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    ht = abs(h.xdim*h.xsize);
    wt = abs(h.ydim*h.ysize);
    biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.01    0.05    wt*scl    ht*scl]);
    
    fig2=subplot (222);
    imagesc(squeeze(d(:,y,:))'), axis ([0.5 h.xdim 0.5 h.zdim]),axis  xy %, axis tight
    caxis([cmin cmax])
    hold on; plot(x,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    %set(gca,'Position',[ 0.5    0.5    0.45    0.45]);
    wt = abs(h.xdim*h.xsize);
    ht = abs(h.zdim*h.zsize);
    %biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.51    0.55    wt*scl    ht*scl]);

    fig3=subplot (221);
    imagesc(squeeze(d(x,:,:))'), axis ([0.5 h.ydim 0.5 h.zdim]), axis xy%, axis tight
    caxis([cmin cmax])
    hold on; plot(y,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
    ht = abs(h.zdim*h.zsize);
    wt = abs(h.ydim*h.ysize);
    %biggest = max([ht wt]);
    scl = 0.45 / biggest;
    set(gca,'Position',[ 0.01    0.55    wt*scl    ht*scl]);
    
end

subplot(224),
%imagesc(0),
caxis([cmin cmax])
%set(gca,'Position',[ 0.5    0.01    0.45  0.45]);
colorbar('West')
axis off
return
