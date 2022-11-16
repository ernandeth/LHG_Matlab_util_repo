
%default values:
%{
mods.kcap = 10; %2;
mods.kgl = 10;
mods.karf = 50;
mods.kel = 50;
mods.kgna = 5; %10;
mods.kgk = 10;
%}

figure(5)
USflag = 2;  % G_leakage modulation
USflag = 5;  % G_Na modulation
USflag = 1;  % Cap. modulation

prf_vals = [0.01:0.05:1 1:10 ];
APmax_prf = [];

for n=1:length(prf_vals)
    
    prf=prf_vals(n);
    
    HHrk4_prf
    drawnow
    
    APmax_prf = [APmax_prf max(vsave(:,1))]
end

figure(105); plot(prf_vals, APmax_prf), title('Modulation of G\_C_m with different PRF')
xlabel('PRF (kHz)')
ylabel('Peak membrane Voltage')
print mod_prf_cap -dpng

