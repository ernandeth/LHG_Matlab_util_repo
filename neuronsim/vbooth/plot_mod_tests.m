myDirs = dir('cf_*')

peaks=[];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of gl modulation
    load test_mods.mat
    peaks=[peaks; APmax_kgl];

    cd ..
end
plot(kgl_vals, peaks)
title('modulation of Leakage Conductance')
xlabel('modulation factor (kgl)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')

%%

figure
peaks=[];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of gl modulation
    load test_mods.mat
    peaks = [peaks; APmax_kgna];

    cd ..
end
plot(kgna_vals, peaks)
title('modulation of Sodium Conductance')
xlabel('modulation factor (kgl)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')

%%

figure
peaks = [];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of Capacitance modulation
    load test_mods.mat
    peaks = [peaks; APmax_kcap];
    
    cd ..
end

plot(kcap_vals, peaks);
title('modulation of Capacitance')
xlabel('modulation factor (kgl)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')

%%

figure
peaks = [];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of Capacitance modulation
    load test_mods.mat
    peaks = [peaks; APmax_kgk];
    cd ..
end
plot(kgk_vals, peaks);
title('modulation of Potassium Conductance')
xlabel('modulation factor (kgk)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')


%%
figure
peaks = [];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of modulation of the ARF
    load test_mods.mat
    peaks = [peaks; APmax_karf];
    cd ..
end
plot(karf_vals, peaks);
    
title('modulation of Acoustic Radiation Force')
xlabel('modulation factor (karf)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')

%%

figure
peaks = [];
for d=1:length(myDirs)
    cd(myDirs(d).name)
    % making a plot of El modulation
    load test_mods.mat
    peaks = [peaks; APmax_kel];
    cd ..
end

plot(kel_vals,peaks);
title('modulation of Reversal Potential')
xlabel('modulation factor (kel)')
legend( '100 KHz', '500 KHz','1000 KHz', '2500 KHz', '5000 KHz')
