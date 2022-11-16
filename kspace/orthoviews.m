function orthoviews(im,dolog,bxlabel,bylabel,bzlabel)

    if nargin<5 || isempty(bzlabel)
        bzlabel = ["" ""];
    end
    
    if nargin<4 || isempty(bylabel)
        bylabel = ["" ""];
    end
    
    if nargin<3 || isempty(bxlabel)
        bxlabel = ["" ""];
    end
    
    if nargin<2 || isempty(dolog)
        dolog = 0;
    end

    if nargin<1 || isempty(im)
        figure, im = lightbox('timeseries_mag'); close gcf
    end
    
    if sum(imag(im)>0,'all'), im = abs(im); fprintf('iscomplex'), end
    if dolog, im = 1e3*log(im+eps()); end
    im = (im - min(im,[],'all'))/(max(im,[],'all') - min(im,[],'all'));
    
    % Plot cut through center of first dimension
    subplot(1,3,1)
        imshow(squeeze(abs(im(round(end/2),end:-1:1,end:-1:1)))')
        axis on
        xlabel('Y'); ylabel('Z');
        xticks([1 size(im,2)]), xtickangle(0)
        xticklabels(bylabel)
        yticks([1 size(im,3)]), ytickangle(90)
        yticklabels(bzlabel(end:-1:1))
        
    % Plot cut through center of second dimension
    subplot(1,3,2)
        imshow(squeeze(abs(im(end:-1:1,round(end/2),end:-1:1)))')
        axis on
        xlabel('X'); ylabel('Z');
        xticks([1 size(im,1)]), xtickangle(0)
        xticklabels(bxlabel)
        yticks([1 size(im,3)]), ytickangle(90)
        yticklabels(bzlabel(end:-1:1))
        
    % Plot cut through center of third dimension
    subplot(1,3,3)
        imshow(squeeze(abs(im(end:-1:1,end:-1:1,round(end/2))))')
        axis on
        xlabel('X'); ylabel('Y');
        xticks([1 size(im,1)]), xtickangle(0)
        xticklabels(bxlabel)
        yticks([1 size(im,2)]), ytickangle(90)
        yticklabels(bylabel(end:-1:1))
    
end

