clear; clc;

fids = dir('chirp*.fid');
gamma = 4257.748;  % Hz / Gauss




allG = [];
allInputs = [];
allrealp = [];
allimagp = [];
allH_nice = [];
allH2_nice = [];
allgirf = [];

for n=5:7
    % 5: coronal -> z
    % 6: axial  -> y
    % 7: sagittal -> x
    
    % Import and process data
    [ A procpar] = getFid(fids(n).name);
    
    % only the first 10 K points are good... too noisy afterward
%     NPOINTS = 5e3;
%     A = A(1:NPOINTS,:);
    
    [a b] = size(A);
        
    % every other echo has a gradient trajectory turned on
    % first half of signals = reference FID
    % second half = FID in presence of gradient
    idx1 = 1:b/2;
    idx2 = b/2+1:b;
    
    dt = 1/procpar.sw;
    pos = procpar.pss;
    phiref = A(:,idx2);
    phiRaw = A(:,idx1);
    
    % error in z axis: there are only four positions, not six
    pos = pos(1:4);
    
    % Subtract reference phase
    tmp = (phiRaw ./ phiref );
    phi = zeros(size(phiref));
    
    for j=1:4
        phi(:,j) = phase(tmp(:,j));
    end
    
   
    Gis = zeros(length(phi) , length(pos));
    count = 1;
    % for every slice collected 
    for m = 1:length(pos)
        
        % calculate the derivative in time at each time point
        for p = 2:length(phi)-1
            
            Gis(p, count) = ( phi(p+1, m) - phi(p-1, m)) / (gamma*dt*4*pi*pos(m));
        end
        
        Gis(1, count) = (phi(2,m) - phi(1,m)) / (gamma*dt*4*pi*pos(m));
        Gis(end, count) = (phi(end,m) - phi(end-1,m)) / (gamma*dt*4*pi*pos(m));
        count = count + 1;
    end
    allpos = Gis;
    
    NPOINTS = length(Gis);
    
    % let's work with the average Gradient waveform only
    % (average all the slices) we will have ONE gradient waveform per axis
    
    %Gis = Gis(:,4);
    Gis = mean(Gis,2);
    
    %  Insert some mechanism to exclude outlier slices and average the others
    %  together
    
    
    allG = [allG  Gis];
    
    % Get th nominal gradient wavefor:  a 'chirp' function :
    [inputChirp tsamp] = getInputChirp(procpar);
    inputChirp = inputChirp(1:NPOINTS);

    allInputs = [allInputs  inputChirp'];
    
   
end

% There is a 40 ms. delay in the gradient waveform
% throw out the first 40 ms.
% hence, we don't get the last 40 ms of the nominal Chirp either

aqdelay = 40e-6 / dt;
aqdelay = 0;  % or ... may be not.

allG = allG(aqdelay+1:end , :);
allInputs = allInputs(1 : end-aqdelay, :);

% Now get rid of another 100 points because of ring-down from prephasing
% gradient.
allG = allG(4*aqdelay+1 : end , :);
allInputs = allInputs(4*aqdelay+1 : end, :);

% calculate frequency domain of input
allInputsF = fft(allInputs,[],1);

% calculate frequency content of measured output
allGF = fft(allG,[],1);

% Calculate transfer function in frequency domain
allH = allGF ./ (eps + allInputsF);


% a time vector
Npoints = size(allG,1);

t = [0:Npoints-1]*dt;
% Nyquist frequency for the input waveform
% Note that this is not the maximum frequency in the chirp/
Nyq = procpar.sw/2;
df = 2*Nyq/Npoints;
f = linspace(-Nyq, Nyq-df, Npoints);



% we;re only interested in the stuff we played.
% cutoff frequency is the max frequency of the input chirp
cutoff = procpar.f1;
cutoff_N = floor(cutoff * (Npoints/2 -1) / Nyq);

% f2 and H2 are the lower portion of the spectrum: 
% where our input chirp is
f2 = linspace(-cutoff, cutoff-df, cutoff_N*2);
f2pos = f2(end/2:end);

allH2 = [...
    allH(1:cutoff_N,:) ;
    allH(end-cutoff_N+1 : end, :) ];


% for every gradient waveform (they are all stacked slices and axes)
for n=1:size(allG,2)
    
    tmp = allG(:,n);
    
    % Show the original data in time and freq. domains
    figure (1)
    subplot(211)
    plot(t, tmp)
    hold on
    plot(t, allInputs(:,1),'r')
    hold off
    axis tight
    title('Gradient waveforms ')
    legend( 'Measured Output', 'Nominal Input')
    hold off
    
    subplot(212)
    plot(f, fftshift(real(allGF(:,n))))
    hold on
    plot(f, fftshift(real(allInputsF(:,n))),'r')
    hold off
    axis tight
    title('Gradient FFT (real) ')
    legend('Measured Output' , 'Nominal Input')
    
    hold off
     
    drawnow
    
    
    % Now try to fit a polynomial to each transfer function:
    % real and imaginary separately
    % Real part:
    tmp = real(allH2(:,n))  ;
    
    % take the positive part of the spectrum only
    tmp = tmp(1:end/2);
    
    [p,s] = polyfit(linspace(0, cutoff-df ,length(tmp)), tmp', 20);
    allrealp = [allrealp ; p];
    realjunk = polyval(p,linspace(0, cutoff-df ,length(tmp)));
    
    % rebuild the negative part of the spectrum
    % remember:  I havent' done any fftshifts
    realjunk = [ realjunk  0 realjunk(end:-1:2)];

    
    %  Imaginary part:
    tmp = imag(allH2(:,n))  ;   
    tmp = tmp(1:end/2);

    
    [p,s] = polyfit(linspace(0, cutoff-df ,length(tmp)), tmp', 20);
    imjunk = polyval(p,linspace(0,cutoff-df, length(tmp)));
    allimagp = [allimagp ; p];
    
    % rebuild the negative part of the spectrum
    % remember:  I havent' done any fftshifts
    imjunk = [imjunk 0 -imjunk(end:-1:2) ];

    
    % the 'ideal' transfer function:
    H2_nice = complex(realjunk, imjunk);
    
    % stack the responses of all coil combos into a single matrix
    allH2_nice = [allH2_nice ; H2_nice];
    
    tmp = (allH2(:,n))  ;
    
    figure (2)
    subplot(211)
    plot(f2, fftshift(real(tmp)),'k')
    hold on
    plot(f2, fftshift(real(H2_nice)))
    xlabel('Frequency (Hz)')
    hold off
    %axis([0 2000 0 1])
    title('Real of Transfer Function')
    legend('Measured', 'Estimated')
    
    subplot(212)
    plot(f2, fftshift(imag(tmp)),'k');
    hold on
    plot(f2, fftshift(imag(H2_nice)))
    xlabel('Frequency (Hz)')
    hold off
    title('Imag of Transfer Function')
    legend('Measured', 'Estimated')
    
    figure (8)
    subplot(211)
    plot(f2, fftshift(abs(tmp)),'k')
    hold on
    plot(f2, fftshift(abs(H2_nice)))
    xlabel('Frequency (Hz)')
    hold off
    %axis([0 2000 0 1])
    title('mag of Transfer Function')
    legend('Measured', 'Estimated')
    
    subplot(212)
    plot(f2, fftshift(angle(tmp)),'k');
    hold on
    plot(f2, fftshift(angle(H2_nice)))
    xlabel('Frequency (Hz)')
    hold off
    title('Phase of Transfer Function')
    legend('Measured', 'Estimated')
    
    
    % upsample the Transfer function to the original resolution:
    zpad = zeros(1, Npoints-length(H2_nice));
    H_nice = [H2_nice(1:end/2) zpad H2_nice(end/2+1 : end) ];
        
    % the gradient's impulse response function in the imte dpmain
    % is the iFFT of the transfer function in frequency
    girf = ifft(H_nice);
    girf = girf(1:Npoints/2-1);
    
    allgirf = [allgirf ; girf];
    
    % try to predict the output given the input and the impulse response
    % function:
    output_nice = conv(allInputs(:,n)', girf) ;
    output_nice = output_nice(1:Npoints) ;
    
    % Plot outputs and inputs.
    t = linspace(0, procpar.at, Npoints)';
    tmp = allG(:,n);
    
    figure (3)
    subplot(211)
    plot(t, real(output_nice),'k')
    hold on
    plot(t, real(tmp))
    hold off
    axis([0 0.05 -1 1])
    title(['Real of output - axis n. ' num2str(n)] )
    legend('Predicted', 'Measured ')
    
    subplot(212)
    plot(t, imag(output_nice),'k')
    hold on
    plot(t, imag(tmp))
    title('Imaginary of output')
    legend('Predicted', 'Measured ')
    axis([0 0.05 -1e-3 1e-3])
    hold off
    drawnow
    
    
end

% save /Users/hernan/data/girf/girf_parms allrealp allimagp allH_nice dt Nyq allgirf cutoff
