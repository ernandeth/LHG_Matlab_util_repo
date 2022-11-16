
%default values:
mods.kcap = 2;
mods.kgl = 10;
mods.karf = 50;
mods.kel = 50;
mods.kgna = 10;
mods.kgk = 10;

loadResults = 0;

if loadResults
    load test_mods.mat
end

figure
%
%% vary the amount of modulation to the capacitance
if ~loadResults
    figure(1)
    USflag = 1
    kcap_vals = [0.6:0.1:1.4 ];
    APmax_kcap = [];
    
    for kcap =  kcap_vals
        mods.kcap = kcap;
        
        
        HHrk4_prf
        drawnow
        APmax_kcap = [APmax_kcap max(vsave(:,1))];
    end
end

figure(101); plot(kcap_vals, APmax_kcap), title('Modulation of Capacitance');
print mod_cap -dpng
save test_mods_cap.mat ...
    kcap_vals APmax_kcap 
%
figure

%% vary the amount of modulation to the leak conductance
if ~loadResults
    figure(2)
    USflag = 2
    kgl_vals = [2:0.5:6];
    APmax_kgl = [];
    
    for kgl =  kgl_vals
        mods.kgl = kgl;
        
        HHrk4_prf
        drawnow
        
        APmax_kgl = [APmax_kgl max(vsave(:,1))];
    end
end

figure(102); plot(kgl_vals, APmax_kgl), title('Modulation of Leakage conductance')
print mod_gl -dpng
save test_mods_gl.mat ...
    kgl_vals APmax_kgl ...

figure
%{
%% vary the amount of modulation to the ARF
if ~loadResults
    figure(3)
    USflag = 3;
    karf_vals = [0:2:20];
    APmax_karf = [];
    
    for karf =  karf_vals
        mods.karf = karf;
        
        HHrk4_prf
        drawnow
        
        APmax_karf = [APmax_karf max(vsave(:,1))];
    end
end

figure(103); plot(karf_vals, APmax_karf), title('Modulation of ARF')
print mod_arf -dpng
save test_mods_arf.mat ...
    karf_vals APmax_karf

figure


%% vary the amount of modulation to the Reversal Potiential
if ~loadResults
    figure(4)
    USflag = 4;
    kel_vals = [0:2:20];
    APmax_kel = [];
    
    for kel =  kel_vals
        mods.kel = kel;
        
        HHrk4_prf
        drawnow
        
        APmax_kel = [APmax_kel max(vsave(:,1))];
    end
end

figure(104); plot(kel_vals, APmax_kel), title('Modulation of reversal Potential')
print mod_el -dpng
save test_mods_el.mat ...
    kel_vals APmax_kel ...


figure
%}

%% vary the amount of modulation to the Sodium conductance
if ~loadResults
    figure(5)
    USflag = 5;
    kgna_vals = [2.6:0.1:5.2 ];
    %kgna_vals = [20:5:40];
    APmax_kgna = [];
    
    for kgna =  kgna_vals
        mods.kgna = kgna;
        HHrk4_prf
        drawnow
        
        
        APmax_kgna = [APmax_kgna max(vsave(:,1))]
    end
end
figure(105); plot(kgna_vals, APmax_kgna), title('Modulation of Na conductance')
print mod_gna -dpng
save test_mods_gna.mat ...
    kgna_vals APmax_kgna ...

figure
%{
%% vary the amount of modulation to the K conductance
if ~loadResults
    figure(6)
    USflag = 6;
    kgk_vals = [0:2:20];
    APmax_kgk = [];
    
    for kgk =  kgk_vals
        mods.kgk = kgk;
        
        HHrk4_prf
        drawnow
        
        APmax_kgk = [APmax_kgk max(vsave(:,1))];
    end
end
figure(106); plot(kgk_vals, APmax_kgk), title('Modulation of K conductance')
print mod_gk -dpng
save test_mods_gk.mat ...
    kgk_vals APmax_kgk 
%}
save test_mods.mat cf prf ...
    kgl_vals APmax_kgl ...
    kcap_vals APmax_kcap ...
    kgk_vals APmax_kgk ...
    kgna_vals APmax_kgna ...
    kel_vals APmax_kel ...
    karf_vals APmax_karf
%
%

return
    
% making a plot
load cf_0100/test_mods.mat
APmax_kgl1000 = APmax_kgl;
plot(kgl_vals, APmax_kgl1000,'g');
hold on

load cf_0500/test_mods.mat
APmax_kgl700 = APmax_kgl;
plot(kgl_vals, APmax_kgl700,'r');

load cf_1000/test_mods.mat
APmax_kgl600 = APmax_kgl;
plot(kgl_vals, APmax_kgl600,'c');

load cf_2500/test_mods.mat
APmax_kgl500 = APmax_kgl;
plot(kgl_vals, APmax_kgl500,'b');

load cf_5000/test_mods.mat
APmax_kgl100 = APmax_kgl;
plot(kgl_vals, APmax_kgl100,'k');

title('modulation of Leakage Conductance')
xlabel('modulation factor (kgl)')
legend('100 KHz', '500 KHz', '1000 KHz','2500 KHz', '5000 KHz')

