t_seg = 1000;
t_seg_xtra = 1000;
M_PI = pi;
vsasl_venc = 4;
GAMMA_H1 = 26754;

t1 = 0.0;                                               
t2 = (t1 + t_seg) ;
t3 = (t2 + t_seg_xtra + 2*t_seg) ;         
t4 = (t3 + t_seg);
t5 = (t4 + 2*t_seg_xtra + 2*t_seg) ;       
t6 = (t5 + t_seg) ;
t7 = (t6 + t_seg_xtra + 2*t_seg) ;         
t8 = (t7 + t_seg) ;

t1 = t1*1.0e-6;
t2 = t2*1.0e-6;
t3 = t3*1.0e-6;
t4 = t4*1.0e-6;
t5 = t5*1.0e-6;
t6 = t6*1.0e-6;
t7 = t7*1.0e-6;
t8 = t8*1.0e-6;


gvenc_amp = 2.0*2.0*M_PI / (...
    vsasl_venc * GAMMA_H1/(2.0*M_PI) *...
    (t2*t2 - t1*t1 +...
    t4*t4 - t3*t3  +...
    t5*t5 - t6*t6  +...
    t7*t7 - t8*t8)...
    )
