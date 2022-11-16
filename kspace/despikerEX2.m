function indices=despikerEX2(pfilename , showData)
% Usage ... indices = despikerEX(pfile, showData)
%
% looks for spikes in the kspace data and replaces them by the 
% average of neighboring time points
% in all the frames of the indicated Pfile.
% 
% This version finds spikes in the derivatives over time.
%
% this version determines spikes from subtractions of adjacent
% time points 
%
% P file is organized as follows:
%  coil 1:
%  slice 1: baseline interleaves time1, interleaves time 2, etc...
%  slice 2: baseline interleaves time1, interleaves time 2, etc...
%  etc...
%  coil 2: same as above but some random skip between them
%  if #leaves and #points is odd then there is an extra baseline at end
%
%  
%  Slice 1 Frame 0 is baseline (example, expect more)
%  Note: there is frame 0 but not slice 0!
%
%  For sf3d - slc is the phase encode of interest and frm
%  would be the spiral shot of interest (0 is the baseline).
%
%  nslices=info1(3);
%  nechoes=info1(4);
%  nexcitations=info1(5);
%  nframes=info1(6);
%  framesize=info1(9);
%  psize=info1(10);

fastrecflag = 0;

if(nargin<2)
    tmp= abs(fft(fftshift(readraw(pfilename,1)-readraw(pfilename,2))))+...
        abs(fft(fftshift(readraw(pfilename,3)-readraw(pfilename,4))))+...
        abs(fft(fftshift(readraw(pfilename,5)-readraw(pfilename,6))))+...
        abs(fft(fftshift(readraw(pfilename,7)-readraw(pfilename,8))));
    figure,plot(tmp);
    maxpoint=find(tmp==max(tmp(:)))
    bad_range=[maxpoint-2:maxpoint+2]
    drawnow, pause
end


RAWHEADERSIZE=61464;
% old P files :% RAWHEADERSIZE=40080;


HEADERSIZE=61464;
% odl Pfiles:  HEADERSIZE=40080;

RDBRAWHEADEROFFSET=30;  % For version 4 use 30.
LOC1 = RDBRAWHEADEROFFSET + 34;
LOC2 = RDBRAWHEADEROFFSET + 70;
LOC3 = RDBRAWHEADEROFFSET + 186;

if isstr(pfilename),
    fid=fopen(pfilename,'r','l');
    if (fid<3), disp(['Could not open file! ...']); end;
else
    fid=pfilename;
end;

hdr=fread(fid, HEADERSIZE, 'uint8');

status = fseek(fid,LOC1,'bof');
if (status), disp(['Could not seek to file location 1! ...']); end;
info1 = fread(fid,10,'short');

status = fseek(fid,LOC2,'bof');
if (status), disp(['Could not seek to file location 2! ...']); end;
info2 = fread(fid,6,'short');

status = fseek(fid,LOC3,'bof');
if (status), disp(['Could not seek to file location 3! ...']); end;
info3 = fread(fid,50,'float');


nslices = info1(3);
nechoes = info1(4);
nexcitations = info1(5);
nframes = info1(6);
framesize = info1(9);
psize = info1(10);    

nfids = nslices * nechoes * (nframes+1)-1;
nfidsperslice = nechoes * (nframes+1)-1;

status=fseek(fid,RAWHEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;

% data type here
if (psize==4), 
    fmt=['long']; 
else, 
    fmt=['short']; 
end;

indices=[];

% open the output file and stick the header and baseline in it
ofp=fopen(sprintf('f_%s', pfilename) ,'wb', 'ieee-le');
fwrite(ofp, hdr, 'uint8');

% do one slice at a time
for sl=1:nslices
    
    % read the complex data in  
    baseline=fread(fid, 2*framesize, fmt);

    raw = fread(fid, [2*framesize nfidsperslice], fmt);
    raw = raw';
    re_raw = raw(:, 1:2:2*framesize);
    im_raw = raw(:, 2:2:2*framesize);
    raw = complex(re_raw, im_raw);
    
    
    %end;
    

    
    if nargin<2, showData=0; end;

    % previously - if fastereceiver was on, we had to demodulate the data
     
     if (fastrecflag),
         %demodulation?
         th = 2*pi*[1:framesize]'/4;
         th = kron(th , ones(1,size(raw,1)));
         draw = raw.*exp(-i*th');
     else
         draw =raw;
     end;
    

    out = zeros(size(draw));
    out=draw;
    % now do the filter (very crude stuff)

    means = mean(out,1);
    stds2 = std(out,0,1);
 
    % filter 2:  derivatives
     derivs = (diff(abs(out),1,1));
     stds = std(derivs,0,1);

     % process the time series of each point on the FID
     for col=1:size(out,2)
         badrows = find((derivs(:,col) > 2*stds(col)) );
         fprintf('\rProcessing point %d of %d slice %d',col,size(out,2),sl);
         for r=1:length(badrows)
             row=badrows(r);
             if (abs(out(row,col))> abs(out(row+1,col) ))
                 out(row,col)=means(col);
             else
                 out(row+1,col)=means(col);
            end
         end
     end
     
     
    
    if (showData==1),
        %if (nargout==0),
        if (nargin==1),
            info1
            info2
            info3
        else,
            subplot(221)
            imagesc(abs(raw));
            title('input')
            colorbar
            subplot(222)
            imagesc(abs(out));
            title('output')
            colorbar
            subplot(223)
            imagesc(abs(draw - out));
            title('abs(input - output)')
            colorbar
            for ct=2:size(raw,1)
                subplot(224)
                %plot(out2(ct,1:2:end));
                %title(sprintf('output FID number: %d',ct));
                %hold on
                %plot(out2(ct,2:2:end),'r');
                %hold off
                plot(abs(draw(ct-1,:) - draw(ct,:)));
                hold on
                plot(abs(out(ct-1,:) - out(ct,:)),'r');
                hold off
                title('(thisFID - lastFID)'); 
                legend('in', 'out')
                %axis([0 4000 0 100])
                
                %pause
                fprintf('FID number: %d ', ct)
                subplot(223)
                plot(abs(draw(ct,:)));
                hold on
                plot(abs(out(ct,:)),'r');
                hold off
                title(sprintf('input FID # %d',ct));
                legend('in', 'out')
                %axis([0 4000 0 100])
                drawnow
            end;
        end
    end
    
     if (fastrecflag),
         % modulation?
         out = out.*exp(i*th');
     end;
    %%%%%%%%%%
    %out=raw;
    %%%%%%%%%%
    
    re_out = real(out);
    im_out = imag(out);
    
    % put the data in the right format for writing
    out2=zeros(nfidsperslice , framesize*2);
    for count=1:framesize
        out2(:,2*count-1) = re_out(:,count);
        out2(:,2*count) = im_out(:,count);
    end
    
    fwrite(ofp, baseline, fmt);
    for ct=1:size(out2,1)
        fprintf('\rWriting point %d of %d ..slice %d',ct,size(out2,2),sl);
        
        fwrite(ofp, out2(ct,:), fmt);
    end
    
end

if isstr(pfilename),
    fclose(fid);
end;
fclose(ofp);


