% Make a simple inversino recovery experiment.
t_del = [0.2 0.2 0.2:0.2:2]'
t_adj = 10*ones(size(t_del));
doAS = zeros(size(t_del));
isSel = doAS;

save t_delays.txt t_del -ascii
save t_adjusts.txt t_adj -ascii
save isVelocitySelective.txt isSel -ascii
save doArtSuppression.txt doAS -ascii

