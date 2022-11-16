% using nominal schedule:
t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');
if exist(fullfile(pwd,'t_aqs.txt'))
	t_aq = load('t_aqs.txt');
else
    t_aq = 0.0349*ones(size(t_tag));
end

TR = t_adjust + t_tag + t_delay + t_aq
hold off
plot(TR)
hold on
plot(t_tag,'r')
plot(t_delay,'k')
plot(t_adjust,'g')
legend('TR', 't tag', 't delay', 't adjust')
