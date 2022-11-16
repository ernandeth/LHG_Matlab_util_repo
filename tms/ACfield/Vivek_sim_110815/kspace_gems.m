function [kro,kpe,xx_mat,yy_mat] =kspace_gems(procpar)


%--------------- Calculation for kro-------------------------------------%

kro_line = zeros(1,(procpar.np / 2)); % Memory allocation for each k-space line
del_kro = procpar.gambar * procpar.gro * (1 / procpar.sw) ; % Step in k-space = 1/FOV_ro

if strcmp(procpar.orient,'trans')
%     kro_start = - (procpar.gambar * -procpar.gror * procpar.tpe); % Using value of dephasing gradient
        kro_start = - ( procpar.gambar * -procpar.gro * procpar.at) / 2; % Using half of readout
    %     gradient
    
    for timepoint = 1:1:(procpar.np/2) % Calculation for a single line in k-space
        if timepoint  == 1
            kro_line(timepoint) = kro_start;
        else
            
            kro_line(timepoint) = kro_line(timepoint - 1) - del_kro;
        end
    end
    
end

if strcmp(procpar.orient,'cor')
%     kro_start = (procpar.gambar * procpar.gror * procpar.tpe); % Using value of dephasing gradient
        kro_start = ( procpar.gambar * procpar.gro * procpar.at) / 2; % Using half of readout gradient
    
    for timepoint = 1:1:(procpar.np/2) % Calculation for a single line in k-space
        if timepoint  == 1
            kro_line(timepoint) = kro_start;
        else
            
            kro_line(timepoint) = kro_line(timepoint - 1) - del_kro;
        end
    end
    
end




if strcmp(procpar.orient,'sag')
%     kro_start = (procpar.gambar * procpar.gror * procpar.tpe); % Using value of dephasing gradient
        kro_start = ( procpar.gambar * procpar.gro * procpar.at) / 2; % Using half of readout gradient
    
    for timepoint = 1:1:(procpar.np/2) % Calculation for a single line in k-space
        if timepoint  == 1
            kro_line(timepoint) = kro_start;
        else
            
            kro_line(timepoint) = kro_line(timepoint - 1) - del_kro;
        end
    end
    
end


kro = repmat(kro_line,1,(procpar.np/2)); % Each readout line in k-space is repeated procpar.nv times

%-------------------------------------------------------------------------%



%--------------------------Caculation for kpe------------------------------

kpe = zeros(procpar.np/2,procpar.nv);
% del_kpe = 1/procpar.lpe;
del_gpe = 1 / ( procpar.gambar * procpar.tpe * procpar.lpe );

% pe_gradvalue = del_gpe .* -( (1:procpar.nv) - (procpar.nv+1)/2);
% 
pe_gradvalue = del_gpe .* -procpar.petable;

kpe_line = -(procpar.gambar .* pe_gradvalue * procpar.tpe)'; % Kpe(t) is the area under phase encode gradient

for pestep = 1:procpar.nv
    kpe(pestep,:) = kpe_line(pestep);
end
kpe = kpe';
kpe = kpe(:)';


del_x = procpar.lpe / procpar.nv;
del_y = procpar.lro / (procpar.np/2);

x = linspace(-(procpar.nv/2),(procpar.nv/2)-1,procpar.nv)*del_x;
y = linspace(-(procpar.np/4),(procpar.np/4)-1,(procpar.np/2))*del_y;
[xx_mat yy_mat] = meshgrid(x,y);
yy_mat = -yy_mat;
xx_mat = -xx_mat;

return

% %--------------------------------------------------------------------------

% Test K-space for axial sections

start = tic;

n = procpar.nv;
orig_signal= ReadVarian2D(filename_fid);
orig_signal = reshape(orig_signal,n,n);
orig_signal  = fliplr(flipud(orig_signal));
m = ifftshift(ifft2(fftshift(orig_signal)));
% m=phantom(n,n);
signal_estimate = zeros(1,n*n);


del_x = procpar.lpe / procpar.nv;
del_y = procpar.lro / (procpar.np/2);

x = linspace(-(n/2),(n/2)-1,n)*del_x;
y = linspace(-(n/2),(n/2)-1,n)*del_y;
[xx_mat yy_mat] = meshgrid(x,y);
yy_mat = -yy_mat; % Vnmrj axes
xx_mat = -xx_mat;
% xx = xx_mat(:);
% yy = yy_mat(:);

for itime = 1: n*n
    
%         -----------------------------------Calculation for signal-------------
    
    signal_estimate(itime) = sum(sum(m .* exp(-1i * 2 * pi *(yy_mat .* kro(itime)' + xx_mat .* kpe(itime)))));
    
end

% toc(start);

signal = reshape(signal_estimate,n,n);
% signal = fliplr(flipud(signal));
image1 = ifftshift(ifft2(fftshift(signal)));
figure(),subplot(321),imagesc(abs(image1)),colormap(gray),
subplot(322),imagesc(abs(m)),colormap(gray)
% subplot(322),imagesc(abs(ifftshift(ifft2(fftshift(orig_signal))))),colormap(gray)
subplot(323),imagesc(abs(signal)),colormap(gray)
subplot(324),imagesc(abs(orig_signal)),colormap(gray)
subplot(325),imagesc(angle(signal)),colormap(gray)
subplot(326),imagesc(angle(orig_signal)),colormap(gray);

figure;plot(1:length(signal(:)),phase(signal(:)),'r',... 
     1:length(signal(:)),phase(orig_signal(:)),'g');

% count_row = 0;
% count_col = 0;
% 
% 
% signal_shift = zeros(procpar.nv,procpar.np/2);
% temp = zeros(procpar.nv,procpar.np/2);
% for ii = 1:procpar.nv
%     row_increment = ii+1;
%     if row_increment > procpar.nv
%         row_increment = 1+count_row;
%         count_row=count_row+1;
%     end
%     
%     
%     temp(row_increment,:) = signal(ii,:);
%     
% end
% 
% for jj = 1:procpar.np/2
%     col_increment = jj+1;
%     if col_increment > procpar.np/2
%         col_increment = 1+count_col;
%         count_col=count_col+1;
%     end
%     signal_shift(:,col_increment) = temp(:,jj);
%     
% end
%     
% 
% imagesc(abs(signal_shift)),colormap(gray),figure(),imagesc(abs(signal)),colormap(gray),figure(),imagesc(abs(orig_signal)),colormap(gray)
% 
% image1 = ifftshift(ifft2(fftshift(signal_shift)));
% figure(),imagesc(abs(image1)),colormap(gray)







