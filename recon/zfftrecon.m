function zfftrecon(root, nslices, doCenterOut)
% function zfftrecon(root, nslices, doCenterOut)
%
% this program is written for the purpose of doing the
% Z direction FFT in 3D spiral (stack of spirals acq)
% it treats every .img file as a singles slice

% read the phase and magnutide images separetlye and
% make complex numbers out of them:
phs = read_img_series(['p_' root]);
[mag, h] = read_img_series(root);
data = mag.*exp(-i.* phs/1000);


% fix the headers to reflect the size of the whole volume
h.zdim = nslices;
h.tdim = size(data,1)/nslices;

% allocate space for the output
outdata = zeros(h.tdim, h.xdim*h.ydim*h.zdim);
order=1:nslices;
if doCenterOut
    order(1:2:end) = nslices/2 + 1 : nslices
    order(2:2:end) = nslices/2: -1: 1
end

for t=0:h.tdim-1
%keyboard
	% put together all the slices into matrices whose rows
	% represent the z direction:
	buf = data(t*nslices + 1 : (t+1)*nslices,  : );
    if doCenterOut
        tmp = zeros(size(buf));
        for c=1:nslices
            tmp(order(c),:) = buf(c,:);
        end
        buf = tmp;
    end
	% noe do the FFT along the Z direction
	BUF = ifftshift(ifft(ifftshift(buf,1) ,[],1),1).';
	% put them into the output matrix
	outdata(t+1,:) = BUF(:);
end
warning off
% write out the data
write_img('3Dout.img', abs(outdata),h);
write_hdr('3Dout.hdr',h);

% write out the data
write_img('p_3Dout.img', 1000*angle(outdata),h);
write_hdr('p_3Dout.hdr',h);
warning on
return

