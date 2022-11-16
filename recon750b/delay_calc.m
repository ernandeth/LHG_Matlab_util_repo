function [delayIn,delayOut] = delay_calc(pfile,strShowPlots,sliceToShow)
% delay_calc.m - given Pfile, autodetects delays 
% 
% INPUTS
% pfile - string, name of pfile
% 'plots' - [optional], generates plots showing a slice with different delay terms.
% sliceToShow - if 'plots' specified, have it plot specific slice (defaults to slice 5)
% 
% OUTPUTS
% delayIn - calculated delay for reverse spiral
% delayOut - calculated delay for forward spiral
% 
% EXAMPLES
% [1] [delayIn,delayOut] = delay_calc(pfile)  % Calculates delay terms
% 
% [2] [delayIn,delayOut] = delay_calc(pfile,'plots')  % Calculates delays and shows images for slice #5
% 
% [3] [delayIn,delayOut] = delay_calc(pfile,'plots',8)  % Calculates delays and shows images for slice #8

% $Id: delay_calc.m 675 2013-06-03 16:34:50Z klitinas $

% Default to not showing plots
if ~exist('strShowPlots','var')
    strShowPlots = 'noplot';
end

% Default slice to plot = 5
if ~exist('sliceToShow','var')
    sliceToShow = 5;
end


posdelmax=20;
negdelmax=-20;

sprec1(pfile,'cp','m');

dellist=negdelmax:posdelmax;

hdr = getpfilehdr(pfile);
revflg = hdr.rdb.user21;

% Spiral in
if revflg > 0
    for x=1:length(dellist)
        sprec1(pfile,'cp','h','t','1','QI','d',int2str(dellist(x)));
        tmp=dir('vol*.img');
        fname=tmp.name;
        [img,hdr] = read_img(fname);
        imgR = reshape(img,hdr.xdim,hdr.ydim,hdr.zdim);
        %ddat(:,:,:,x)=double(analyze75read(fname));
        ddat(:,:,:,x)=imgR;
        !rm vol*.img vol*.hdr
    end
    
    % Do plots
    if strcmpi(strShowPlots,'plots')
        strTitle = sprintf('%s, slice %d: spiral in',pfile,sliceToShow);
        inSlice = squeeze(ddat(:,:,sliceToShow,:));
        subplotimg(inSlice,dellist,strTitle)
    end
end
%%%%%%%%%%%%%%%%%%%

% Spiral out
if revflg ~= 1
    for x=1:length(dellist)
        sprec1(pfile,'cp','h','t','1','QO','d',int2str(dellist(x)));
        %sprec2(pfile,'cp','h','Q','t','1','d',int2str(dellist(x)));
        tmp=dir('vol*.img');
        fname=tmp.name;
        [img,hdr] = read_img(fname);
        imgR = reshape(img,hdr.xdim,hdr.ydim,hdr.zdim);
        % rdat(:,:,:,x)=double(analyze75read(fname));
        rdat(:,:,:,x)=imgR;
        !rm vol*.img vol*.hdr
    end
    
    if strcmpi(strShowPlots,'plots')
        outSlice = squeeze(rdat(:,:,sliceToShow,:));
        strTitle = sprintf('%s, slice %d: spiral out',pfile,sliceToShow);
        subplotimg(outSlice,dellist,strTitle)
    end
    
end

if revflg == 2
    for x=1:size(ddat,4)
        for y=1:size(rdat,4)
            temp=(ddat(:,:,:,x)-rdat(:,:,:,y)).^2;
            zdiff(x,y)=sum(temp(:));
            dldiff(x,y)=abs(dellist(x)+dellist(y));
        end;
    end;
    ind=find(dldiff(:)<2);
    [del_in,del_out]=find(zdiff==min(zdiff(ind)));
    
    delayIn = dellist(del_in);
    delayOut = dellist(del_out);
    
    fprintf('\nDelay  In: %d',delayIn);
    fprintf('\nDelay Out: %d\n\n',delayOut);
    
else
    delayIn = [];
    delayOut = [];
end
% inclusive ?
% [m,i] = min(zdiff(:));
% [ix,iy] = ind2sub([41 41],i);
% del_in2 = dellist(ix)
% del_out2 = dellist(iy)
% %