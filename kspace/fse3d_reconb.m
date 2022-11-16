function fse3d_recon2(Pfile , x,y,nechos,nframes)
%function fse3d_recon2(Pfile , xdim, ydim, nechos, nframes)
%
% reconstruction of stack of spirals data using gsp20a for 
% x and y directions and doign the third direction by hand.
% assumes that the P file header info thinks that there is only 
% one slice.
%
str=sprintf('!gsp20a -p %s', Pfile)
eval(str);

outnum=1;
echo=0;
tmp=[];
kspace=[];

for frame = 1:nframes*(nechos)
    p_name=sprintf('sl1.phs.%03d',frame);    
    fprintf('\rreading ... %s  ', p_name);
 
    fp=fopen(p_name,'r');
    p = fread(fp,x*y,'short'); 
    p=p/1000;
    fclose(fp);
            
    m_name=sprintf('sl1.%03d',frame);    
    fprintf('\rreading ... %s  ', m_name);
    
    fp=fopen(m_name,'r');
    m = fread(fp,x*y,'short'); 
    fclose(fp);
    
    tmp = m.*exp(-i.*p);
    % correct the phase cycling here:
    %tmp = (-1)^echo * tmp;
    kspace = [kspace ; tmp'];

    echo=echo +1;

    if echo>=nechos 
    	z = nechos;
    	obj =100* fftshift(ifft(fftshift(kspace,1), [], 1),1);
        phi=zeros(size(obj))l
        for sl=1:nechos
            phi(:,:,sl) = angle(obj(32,32,sl));
        end

    	str=sprintf('!avwcreatehd %d %d %d 1 1.0 1.0 3.0 1.0 %d %d %d 4 vol_%04d.hdr',...
    	    x, y, z, x/2, y/2, z/2, outnum);
    	eval(str);
    	outh = read_hdr(sprintf('vol_%04d.hdr',outnum));
    	fprintf('\r  writing ... vol_%04d.img',outnum);
    	write_img(sprintf('vol_%04d.img', outnum),abs(obj)',outh);
    	
    	str=sprintf('!avwcreatehd %d %d %d 1 1.0 1.0 3.0 1.0 %d %d %d 4 pvol_%04d.hdr',...
    	    x, y, z, x/2, y/2, z/2, outnum);
    	eval(str);
    	fprintf('\r  writing ... pvol_%04d.img',outnum);
    	write_img(sprintf('pvol_%04d.img', outnum),1000*angle(obj)',outh);

    	outnum=outnum+1;
	
	echo =0;
	kspace=[];
    end

end

%!rm sl1.*
return
% renumber the time series....
vols = dir('vol_*.img');
for count=1:length(vols)
    curr_name=vols(count).name;
    curr_name = curr_name(1:end-4);
    str=sprintf('!mvimg %s volume_%04d', curr_name, count)
    eval(str);
    str=sprintf('!mvimg p%s pvolume_%04d', curr_name, count)
    eval(str);
end
