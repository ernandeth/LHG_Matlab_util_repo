function fse3d_recon2(Pfile , x,y,nechos,nframes,centricflag)
%function fse3d_recon(Pfile , xdim, ydim, nechos, nframes, centricflag)
%
% reconstruction of stack of spirals data using gsp20a for 
% x and y directions and doign the third direction by hand.
% assumes that the P file header info thinks that there is only 
% one slice.
%
str=sprintf('!gsp21a -p %s', Pfile)
eval(str);

p_files = dir('*.phs.*');
outnum=1;
for frame = 1:nechos:nframes*(nechos)
    
    outim=[];
    
    if(centricflag)
        encvals=zeros(1,nechos);
        encvals(1:2:end)=nechos/2+1:nechos;
        encvals(2:2:end)=nechos/2:-1:1;
    else
        encvals=1:nechos;
    end

    
    outim=zeros(nechos,x*y);
    
    for echo=0:nechos-1
        fprintf('\rreading ... %s  ', p_files(frame+echo).name);
 
        fp=fopen(p_files(frame+echo).name,'r');
        p = fread(fp,x*y,'short'); 
        fclose(fp);
                
        m_name = p_files(frame+echo).name;
        m_name = [m_name(1:3) m_name(8:end)];
        
        
        fp=fopen(m_name,'r');
        m = fread(fp,x*y,'short'); 
        fclose(fp);
        
        tmp = m.*exp(i*p/1000);
        
        %outim = [outim ; tmp'];
        outim(encvals(echo+1),:)=tmp';
        
%        tp1=x*y*(encvals(echo+1)-1)+1
%        outim(,tp1:tp1+15) = tmp';
       
    end
    
    z = nechos;
    str=sprintf('!avwcreatehd %d %d %d 1 1.0 1.0 3.0 1.0 %d %d %d 4 vol_%04d.hdr',...
        x, y, z, x/2, y/2, z/2, outnum);
    eval(str);
    tmp2 = abs(fftshift(ifft(fftshift(outim,1), [], 1),1));
    %outim = reshape(outim, h.xdim*h.ydim*h.zdim*nechos, 1);
    outh = read_hdr(sprintf('vol_%04d.hdr',outnum));
    fprintf('\r  writing ... vol_%04d.img',outnum);
    write_img(sprintf('vol_%04d.img', outnum),tmp2',outh);
    
    str=sprintf('!avwcreatehd %d %d %d 1 1.0 1.0 3.0 1.0 %d %d %d 4 p_vol_%04d.hdr',...
        x, y, z, x/2, y/2, z/2, outnum);
    eval(str);
    tmp2 = 1000*angle(fftshift(ifft(fftshift(outim,1), [], 1),1));
    %outim = reshape(outim, h.xdim*h.ydim*h.zdim*nechos, 1);
    outh = read_hdr(sprintf('p_vol_%04d.hdr',outnum));
    fprintf('\r  writing ... p_vol_%04d.img',outnum);
    write_img(sprintf('p_vol_%04d.img', outnum),tmp2',outh);
    outnum=outnum+1;

end

!rm sl1.*
return
% renumber the time series....
vols = dir('vol_*.img');
for count=1:length(vols)
    curr_name=vols(count).name;
    curr_name = curr_name(1:end-4);
    str=sprintf('!mvimg %s volume_%04d', curr_name, count)
    eval(str);
    str=sprintf('!mvimg p%s p_volume_%04d', curr_name, count)
    eval(str);
end
