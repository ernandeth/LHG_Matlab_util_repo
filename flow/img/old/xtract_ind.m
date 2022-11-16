
function tcourse = xtract_ind(hdr, data, x, y, z)
%
% function tcourse = xtract_ind(hdr, data, x, y, z)
%
% here the data is a 2D matrix with all the time points in rows
% x,y,z are vectors defining a cubic ROI
%
% the output is another 2D matrix with the individual time courses in
% selected boundaries
%
    num = 0;

	startx = min(x);
	starty = min(y);
	startz = min(z);
    endx = max(x);
	endy = max(y);
    endz = max(z);
    
    if endx>=hdr.xdim, endx=hdr.xdim;  end
    if endy>=hdr.ydim, endy=hdr.ydim;  end
    if endz>=hdr.zdim, endz=hdr.zdim;  end

    if startx <= 1 , startx=1 ; end    
    if starty <= 1 , starty=1 ; end    
    if startz <= 1 , startz=1 ; end    
    
    %whos data
    
    tcourse=[];
    
    for i=startx:endx
          for j=starty:endy
                for k=startz:endz
                    tcourse = [tcourse  ...
                        data( :,(k-1)*(hdr.xdim*hdr.ydim) + (j-1)*hdr.xdim + i-1)];

                end
          end
    end

    fprintf('\n box boundaries...')
    fprintf('from %d, %d, %d to %d, %d, %d',startx, starty, startz, endx, endy,  endz)
    fprintf('... %d voxels in cube ', num)

    
return