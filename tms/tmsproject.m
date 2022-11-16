function result = project(tms_points)

% function result = project(tms_points, volume)
%
% Luis hernandez
% last edit : 8-13-98
%
% projects a set of points onto a volume's surface (user is prompted for
% an Analyzes file).
% the voulme must be a 3D array where all points outside the
% voulme have zero intensity
% The points are a 3 column data set containing xyz coordinates:
% a center point and two other points co-planar to it.
% the normal to that plane is computed and the point is projected along
% that normal.

	% Select and Read Analyze file.

	[n p] = uigetfile('*.img','Select image file with projection target volume');
	imgname = strcat(p,n);
	sz = size(imgname);
	hdrname = strcat(imgname(1,1:sz(2)-4), '.hdr');
	%
	hdr = read_hdr(hdrname)
	slices = read_img(hdr, imgname);
	
	     
     
	cla;
	set(gca,'Color',[0 0 0]);
	plot3(tms_points(:,1), tms_points(:,2), tms_points(:,3), 'bo');	
    

	sz = size(tms_points);
	for i = 1 : 4 : sz(1)

		% Extract three points at a time.  (assume this configuration:)
		%
		%		4
		% 	2		3
		%		1
		%
		% Points 2,3,4 are extracted. point 1 is not used.
		% Point #4 will be the one projected.

		v = tms_points( i+1 : i+3, :);

               
		plot3(v(1,1), v(1,2), v(1,3), 'y*');	
		plot3(v(2,1), v(2,2), v(2,3), 'y*');	
		plot3(v(3,1), v(3,2), v(3,3), 'g*');	
          
          	center = v(3,1:3);
          	vector1 = v(1,1:3) - v(3,1:3);
          	vector2 = v(2,1:3) - v(3,1:3);
          

		% Calculate the UNIT normal vector to the plane

		normal = cross(vector1, vector2);
		magn = sqrt( sum(power(normal,2) ) );
		normal = normal/magn;

		found_surface = 0;
		project = center;
		step = 0.1;
		k = 0; 
          
                 
          while ( found_surface==0  &  k < 1000 )
	
               project = project - normal*step;
               plot3( project(1),project(2),project(3),'y.')
               %result = [result; project] ;

		% scale the trajectory from cm to pixel units
               xpix = ceil( project(1) * 10 /hdr.xsize );
               ypix = ceil( project(2) * 10 /hdr.ysize );    
               zpix = ceil( project(3) * 10 /hdr.zsize );
               
                             
               % make sure that the TMS point is not outside the image volume
               if ( xpix > 0 & xpix <= hdr.xdim )& ...
                  ( ypix > 0 & ypix <= hdr.ydim )& ...
                  ( zpix > 0 & zpix <= hdr.zdim )
               
			if slices( zpix ).data( ypix,xpix ) ~= 0
                	% note that the image data matrix is transposed !
                      	  found_surface = 1;
                         result = [result; project] ;
                            
                            
                         % Check to see if the user wants to expand to 1 cm areas
                         % This is not done when he wants to compute the surface area
                         % and the center of mass.
                         add_neighbors = get(findobj('Tag','ExpansionCheckBox'),'Value');
                            
                         if add_neighbors ==1
                            disp('adding neighbors...');

                            % also add five mm in every direction
                            
                            hdr 
                            
                            x_expand = ceil(5 / hdr.xsize);
                            y_expand = ceil(5 / hdr.ysize);
                            z_expand = ceil(5 / hdr.zsize) ;                            
                            
                             
                            for nx =-x_expand:x_expand
                               for ny=-y_expand:y_expand
                                  for nz=-z_expand:z_expand
                                     
                                     xx = xpix + nx;
                                     yy = ypix + ny;
                                     zz = zpix + nz;
                                     
                                     if slices( zz ).data( yy, xx ) ~= 0
                                           neighbor = [xx * hdr.xsize/10  yy* hdr.ysize/10  zz*hdr.zsize/10];
                                           result = [result; neighbor];
                                     end
                                        
                                  end
                               end
                            end
                            
                            
                         end
                            
                      	 plot3( project(1),project(2),project(3),'ro');
                         % pause      
		end
                 
              end
              
              k = k+1;
              
          end
          
%plot3(v(1,1), v(1,2), v(1,3), 'b*');	
%plot3(v(2,1), v(2,2), v(2,3), 'b*');	
%plot3(v(3,1), v(3,2), v(3,3), 'b*');	
          
   end
   

	% reattach intensity 
	sz = size(result)
   	result = [result zeros(sz(1),1)];
   
   	whos result  
	pause(1)
          
return	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hold off
%cla
%colormap(gray)

%imagesc(slices(zpix).data' )
%hold on
%disp  ( slices( zpix ).data( xpix,ypix )  )
%plot( xpix,ypix,'g.');

%plot(256,256,'yo')
%plot(128,128,'bo')
%plot(1,1,'ro')

%pause



%%%%%%%%%%%%%%%%%%%%%%%%
%hold on
%plot( xpix,ypix,'ro');
%pause

