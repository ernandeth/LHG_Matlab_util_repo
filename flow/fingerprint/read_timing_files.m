

labelcontrol =  load('labelcontrol.txt');
t_tag = load('t_tags.txt');
t_delay = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');

if exist(fullfile(pwd,'t_aqs.txt'))
    t_aq = load('t_aqs.txt');
else
    t_aq = 0.0349*ones(size(t_tag));
end

tr = t_tag+t_delay+t_adjust +t_aq;
sum(tr)
plot(tr)


    timing_parms.order = 1;
    timing_parms.Nlabel_group = 1;

    timing_parms.t_tag = t_tag;
    timing_parms.t_delay = t_delay;
    timing_parms.t_adjust = t_adjust;
    timing_parms.isLabel = labelcontrol;
    timing_parms.t_aq = t_aq;
    