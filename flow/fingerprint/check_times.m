t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');
t_aqs = load('t_aqs.txt');
labelcontrol =  load('labelcontrol.txt');



nom_t_tag = load('nominal/t_tags.txt');
nom_t_delay = load('nominal/t_delays.txt');
nom_t_adjust = load('nominal/t_adjusts.txt');
nom_labelcontrol =  load('nominal/labelcontrol.txt');

% the fix!
%
% t_tag =t_tag(2:end);
% t_delay =t_delay(2:end);
% t_adjust =t_adjust(2:end);
% 
% nom_t_tag =nom_t_tag(1:end-1);
% nom_t_delay =nom_t_delay(1:end-1);
% nom_t_adjust =nom_t_adjust(1:end-1);



subplot(311)
plot(t_adjust); hold on
plot(nom_t_adjust, 'r')
hold off
title('t adj')

subplot(312)
plot(t_tag); hold on
plot(nom_t_tag, 'r')
hold off
title('t tag')

subplot(313)
plot(t_delay); hold on
plot(nom_t_delay, 'r')
hold off
title('t delay')


figure
subplot(311)
plot(t_adjust - nom_t_adjust, 'r')
hold off
title('t adj')

subplot(312)
plot(t_tag - nom_t_tag, 'r')
hold off
title('t tag')

subplot(313)
plot(t_delay - nom_t_delay , 'r')
hold off
title('t delay')