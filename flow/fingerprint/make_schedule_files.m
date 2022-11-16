function make_schedule_files(timing_parms, dirName)
% function make_schedule_files(timing_parms, dirName)

t_tag           = timing_parms.t_tag;
t_delay         = timing_parms.t_delay;
t_adjust        = timing_parms.t_adjusts ;
isLabel         = timing_parms.labelcontrol ;
order           = timing_parms.order;
t_aq            = timing_parms.t_aq ;
AS_delay        = timing_parms.ArtSup_delay;
doAS            = timing_parms.doArtSup;
RO_type         = timing_parms.readout_type;
label_type      = timing_parms.label_type;

Nframes = length(t_delay)

figure(1)

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

%%%%%%%%
% make up for scanner error:
%AS_delay = AS_delay-0.1;
%%%%%%%
!mkdir tmpDir
save tmpDir/t_tags.txt      t_tag -ascii
save tmpDir/t_adjusts.txt   t_adjust -ascii
save tmpDir/t_delays.txt    t_delay -ascii
save tmpDir/isVelocitySelective.txt isLabel -ascii
save tmpDir/AS_delays.txt   AS_delay -ascii
save tmpDir/doArtSuppression.txt    doAS -ascii
save tmpDir/label_type.txt  label_type -ascii
save tmpDir/RO_type.txt     RO_type -ascii
save tmpDir/order.txt       order -ascii
save tmpDir/RO_time.txt     t_aq -ascii

str = ['!rm -r ' dirName];
eval(str)
str = ['!mv tmpDir ' dirName]
eval(str)