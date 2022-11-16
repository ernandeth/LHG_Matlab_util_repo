function parms = genStructParmLH_VSI(parm_s)
    parms.mtis0 =     1 ;
    parms.Disp =      40;
    parms.r1blood = 1/1.7;

    parms.f =       parm_s(1);
    parms.cbva =    parm_s(2);
    parms.kfor =    parm_s(3);
    parms.bat2 =    parm_s(4);
    parms.bat =     parm_s(5);
    parms.r1tis =   parm_s(6);
    parms.flip =    parm_s(7);
    
    parms.alpha_ti = parm_s(8);
    parms.alpha_ai = parm_s(9);
    parms.alpha_ts = parm_s(10);
end