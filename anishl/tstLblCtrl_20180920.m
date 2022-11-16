% for control-supercontol hybrid
min = Inf;
for ran2 = 1:1000
    rng(ran2)%4
    Nframes = length(timing_parms_min.t_tag);
    iL2 = (rand(1,Nframes));iL2(iL2<0.33) = -1;iL2(iL2>=0.66) = 1;iL2(iL2<0.66 & iL2>=0.33) = 0;
    timing_parms_min.isLabel = iL2;

    checkCRLB_multi_7parms
    if min>CRLB_opt(3)
        min = CRLB_opt(3);
        iL2min = iL2;
    end
    disp(ran2)
end
    fprintf('the min flow std is: ')
    disp(min)