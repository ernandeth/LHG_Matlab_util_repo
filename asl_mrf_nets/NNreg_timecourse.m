function result = NNreg_timecourse(data, netDir)
% function result = NNreg_timecourse(data, netDir)

data = data / norm(abs(data));

fprintf('\nCalculating perfusion with Neural Net ... ');
str = [netDir '/Mynet_cbf.mat']
load(str)

cbf = predict(Mynet, data) ;
result.cbf = cbf / scale *6000; % adjust units to ml/min/100g

%%
fprintf('\nCalculating Arterial Blood Volume with Neural Net ... ');
str = [netDir '/Mynet_cbv.mat']
load(str)

cbv = predict(Mynet, data) ;
result.cbv = cbv / scale;

%%
fprintf('\nCalculating Bolus Arrival Time with Neural Net ... ');
str = [netDir '/Mynet_bat.mat']
load(str)

bat = predict(Mynet, data) ;
result.bat = bat / scale;

%%
fprintf('\nCalculating T1 relaxation with Neural Net ... ');
str = [netDir '/Mynet_r1.mat']
load(str)

r1 = predict(Mynet, data) ;
result.T1 = 1/r1 * scale;

%%
fprintf('\nCalculating T2 with Neural Net ... ');
str = [netDir '/Mynet_r2.mat']
load(str)

r2 = predict(Mynet, data) ;
result.T2 = 1/r2 * scale;
%%
