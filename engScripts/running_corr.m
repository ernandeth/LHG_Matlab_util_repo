% these are initialized before starting

% ssxx=0;
% ssyy=0;
% ssxy=0;
% ssx2=0;
% ssy2=0;
    

% current data, reference, k is current frame number
	%dat1=tmpdat;
    k=vnum;
    ref1=ref(k);
%
    %tc=dat1;
    tc = asl;
    refe=ref1;
    

                ssxx = ssxx+tc;
                ssx2 = ssx2+tc.^2;
                msxx = ssxx./k;
                msx2 = ssx2./k;
                
                ssyy = ssyy+refe;
                ssy2 = ssy2+refe.^2;
                msyy = ssyy./k;
                msy2 = ssy2./k;
                
                ssxy = ssxy+(refe.*tc);
                msxy = ssxy./k;
                
                sdxx = sqrt(msx2-(msxx).^2);
                sdyy = sqrt(msy2-(msyy).^2);
                
                corrMap=(msxy-(msyy.*msxx))./(sdxx*sdyy);

%                cmaps(:,:,:)=p;
    
%    subplot(2,1,2);
%    show_rt(tile(cmaps));
    %title(strcat('t = ',int2str(k)));
    %colormap(hot);
%drawnow;
%end