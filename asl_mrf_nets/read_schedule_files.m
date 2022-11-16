function timing_parms = read_schedule_files(tmpDir)

tmp = load([ tmpDir '/t_tags.txt']);      
    timing_parms.t_tag = tmp;
tmp = load([ tmpDir '/t_adjusts.txt']);   
timing_parms.t_adjusts =tmp;
tmp = load([ tmpDir '/t_delays.txt']);    
timing_parms.t_delay =tmp;
tmp = load([ tmpDir '/isVelocitySelective.txt']); 
timing_parms.labelcontrol =tmp;
tmp = load([ tmpDir '/AS_delays.txt']);   
timing_parms.ArtSup_delay =tmp;
tmp = load([ tmpDir '/doArtSuppression.txt']);    
timing_parms.doArtSup =tmp;
tmp = load([ tmpDir '/label_type.txt']);  
timing_parms.label_type =char(tmp);
tmp = load([ tmpDir '/RO_type.txt']);     
timing_parms.readout_type =char(tmp);
tmp = load([ tmpDir '/order.txt']);       
timing_parms.order = tmp;
tmp = load([ tmpDir '/RO_time.txt']);     
timing_parms.t_aq = tmp;
if (exist([ tmpDir '/RO_time_scan.txt']))
    tmp = load([ tmpDir '/RO_time_scan.txt']);     
    timing_parms.t_aq = tmp;
end

%% now show the timings in a figure

subplot(321)
plot(timing_parms.t_adjusts); title('pre label delay')
axis tight

subplot(322)
plot(timing_parms.t_delay); title('post label delay')
axis tight

subplot(323)
plot(timing_parms.ArtSup_delay); title('AS delay')
axis tight

subplot(324)
plot(timing_parms.t_tag); title('Label Duration')
axis tight

subplot(325)
stem(timing_parms.labelcontrol); title('Label/Control')
axis tight

subplot(326)
stem(timing_parms.doArtSup); title('Do Art Sup')
axis tight

return
