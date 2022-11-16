
label_dur = sort(repmat( 1.5*rand(40,1) + 0.5 , [1,2])); %  ]TR - 0.600 - 0.300;
PID =       sort(repmat( sin(linspace(0, 2*pi, 40) , [1,2]));
TR = label_dur + PID + 0.3;

timing_parms.PID = PID;
timing_parms.label_dur = label_dur;
timing_parms.TR = TR;