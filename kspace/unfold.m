
close all

[file path]  = uigetfile('*.img','Select magnitude *.img file');  
oldpath = pwd;
cd(path)
sz = size(file);
root = file(1,1:sz(2)-7)

files = dir(strcat(root,'*.img'));
hfiles = dir(strcat(root,'*.hdr'));
sz = size(files);
hfiles(1).name;
hdr = read_hdr(hfiles(1).name);

%UNFOLD simulation
GAMMA = 26752;  % rad/s/Gauss
debug=1;
tseries=[];
NIMGS = size(files,1);
TR = 1;
nyquist = 1/(2*TR);

% the cutoff frequency is expressed in units of Hz.
cutoff = 0.5 * nyquist;

% read all the data into a 2D matrix.
% importantL: these are complex data.
tseries = [];
for n=1:NIMGS
    fprintf('\rReading ...%s',files(n).name);
    mtmp =  read_img_data(hdr, files(n).name);
    
    str=sprintf('p_%s',files(n).name);
    ptmp = read_img_data(hdr, str );
    
    tmp= mtmp .* exp(-i*ptmp/1000);
    tseries=[tseries; tmp];
end

outseries = zeros(size(tseries));
% now apply the filter to all the pixels in the time direction
for n=1:size(tseries,2)
    outseries (:,n)= smoothdata(tseries(:,n), TR, cutoff,5);
end
if (debug)
    smoothdata(tseries(:,200), TR, cutoff,5,1);
end

% now write the series to file
for n=1:NIMGS
    % magnitude 
    str=sprintf('u%s',files(n).name);
    fprintf('\n writing .....%s',str)
    writeim(str,abs(outseries(n,:)) , 'short');
    
    str=sprintf('u%s',hfiles(n).name);
    fprintf('\n writing .....%s',str)
    write_hdr(str,hdr);
    
    
    % phase
    str=sprintf('p_u%s',files(n).name);
    fprintf('\n writing .....%s',str)
    writeim(str,1000*angle(outseries(n,:)) , 'short');
     
    str=sprintf('p_u%s',hfiles(n).name);
    fprintf('\n writing .....%s',str)
    write_hdr(str,hdr);
end

cd(oldpath)
% for n=1:1440
%     str=sprintf('!mv sl1.%03d sl%04d.img',n,n)
%     eval(str);
%     str=sprintf('!cp usl0001.hdr sl%04d.hdr',n)
%     eval(str)
% end
