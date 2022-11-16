function swcall(action)


global raw bw
global spect 
global fid_length ax 
global time_buffer freq_buffer databuffer
global mean_time_buffer mean_freq_buffer

switch (action)
case 'load_display'
    % This is where we load the raw data , and immediately do an FFT
    % all the FIDs are displayed simultaneously as an image
    h = guidata(gcbo);  
    filename = get(h.file_name,'String');
    FIDlength = str2num (get(h.fid_length,'String'));
    num_fids = str2num (get(h.num_fids ,'String'));
    fid_length = str2num (get(h.fid_length,'String'));
    hdr_size = str2num (get(h.hdr_size,'String')) ;  
    bw = str2num (get(h.bw,'String')) ;
        
%    raw = read_cplx(filename, fid_length, num_fids,hdr_size);
    raw = read_cplx2(filename, fid_length, num_fids);
    spect = fftshift(fft(raw));

    time_buffer = raw;
    freq_buffer = spect;
    
    mean_time_buffer = mean(time_buffer,2);
    mean_freq_buffer = mean(freq_buffer,2);
    
    figure
    colormap(gray)
    imagesc(real(raw)');
    
case 'choose_display'
    % This allows the user to chose whether how to display the complex
    % and whether to llok at the time domain data or the Freq. Domain
    h = guidata(gcbo);
    plot_type = get(h.display_menu, 'Value');
    
    
    show_nmr(plot_type,bw);
   
    
       
case 'phase'
    % This allows the user to add phase to the raw data at
    % will and updates the display accordingly.
    
    h = guidata(gcbo);

    phase = 2*pi * get(h.phase, 'Value');

    phased_spect = abs(spect) .* exp( i * (angle(spect) + phase));
    phased_raw = ifft(phased_spect);
    
    time_buffer = phased_raw;
    freq_buffer = phased_spect;
    
    mean_time_buffer = mean(time_buffer,2);
    mean_freq_buffer = mean(freq_buffer,2);
    
    plot_type = get(h.display_menu, 'Value');
    
    ax = axis;
    
    show_nmr(plot_type,bw);
    
    
    
case 'time_course'
    % this extracts a timeseries from the data at
    % the frequency selected by the user.

   h = guidata(gcbo);
   %peak_freq = str2num(get(h.peak_freq , 'String') );
   
   % First show the user a spectrum and let them
   % zoom in and choose a peak
   subplot (211)
   plot_type = get(h.display_menu, 'Value');
   show_nmr(plot_type,bw);
   disp('Zoom in and hit ENTER when finished')

   pause   
      
   coords = ginput(1)
   hold on
   plot(coords(1), coords(2), 'ro');
   hold off
   
   % Exract the time series from teh frequency domain data
   peak_freq = coords(1);
   timeseries = freq_buffer(peak_freq,:);
   subplot (212)
   show_series(plot_type, timeseries);
   subplot(211);
   ax = axis;
   
   
case 'remove_water'
% this function does a fit of the water peak and subtracts the result 
% from the raw data

    h = guidata(gcbo);
    
    init_guess = str2num(get(h.water_guess,'String'));
    fidlength  = size(time_buffer,1);
    
    % use the magnitude of the FID fpr water subtraction
    abs_buffer = abs(time_buffer);
    
    for line=1:size(abs_buffer,2)

        guess =lsqcurvefit('t2_decay', ...
                init_guess ,...
                [1:size(abs_buffer,1)]', ...
                abs_buffer(:,line));
    
    	m0 = guess(1);
	    R2 = guess(2);
	    water = t2_decay([m0 R2],  [1:fidlength]);
        abs_buffer(:,line) = abs_buffer(:,line)- water';
       
   end
    
   freq_buffer =  fftshift(fft(abs_buffer));

   mean_time_buffer = mean(time_buffer,2);
   mean_freq_buffer = mean(freq_buffer,2);
   
   plot_type = get(h.display_menu, 'Value');
   show_nmr(plot_type,bw);
   
   
case 'fit_peaks'
    % Fit the data to clean up the spectrum:
    
    h = guidata(gcbo);
    shifts = str2num(get(h.peak_shifts, 'String'));
    T2s = str2num(get(h.peak_T2s, 'String'));
    
    freq_tmp = freq_buffer;
    time_tmp = freq_buffer;
    
    freq_buffer = fit_spectra(time_buffer);
    time_buffer = ifft(freq_buffer);
    
    plot_type = get(h.display_menu, 'Value');
    show_nmr(7);
  
    hold on
    plot(abs(mean(freq_tmp(fid_length/2:fid_length,:) ,2)), 'r');
    hold off
    whos
    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_nmr(plot_type,bw)
% This plots the average signal in the time series as 
% specified by plot_type
GAMMA = 26752 /(2*pi);  % Hz/Gauss
global time_buffer freq_buffer

 tb = mean(time_buffer,2);
 fb = mean(freq_buffer,2);
 
 ppm_range = bw/128;
 
 fr_axis = linspace(-ppm_range/2, ppm_range/2, size(fb,1)); 
 
 
 switch (plot_type)
    case 1
        plot ( real(tb) );
    case 2
        plot ( imag(tb) );
    case 3
        plot ( abs(tb) );
    case 4
        plot ( angle(tb) );
    case 5
        plot (fr_axis, real(fb) );
    case 6
        plot (fr_axis, imag(fb) );
    case 7
        plot ( fr_axis,abs(fb) );
    case 8
        plot ( fr_axis,angle(fb) );
    end
    axis tight    

return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_series(plot_type,  f_series)
% This plots a time series as specified by plot_type

    switch (plot_type)
 
    case 5
        plot ( real(f_series));
    case 6
        plot ( imag(f_series));
    case 7
        plot ( abs(f_series));
    case 8
        plot ( angle(f_series));
        
    otherwise
    errormesg('Is the Plot Type Correct?');    
end

    axis tight    

return
    
