 
function Qpartials_keplerian(sma, e, inc, ran, ape, m, ps::qLawParams)    
    # Set some constants
    b_petro = ps.b_petro
    m_petro = ps.m_petro
    n_petro = ps.n_petro
    r_petro = ps.r_petro

    # Grab parameters
    sma_t   = ps.oet[1]
    e_t     = ps.oet[2]
    inc_t   = ps.oet[3]
    ran_t   = ps.oet[4]
    ape_t   = ps.oet[5]
    Wsma    = ps.oeW[1]
    We      = ps.oeW[2]
    Winc    = ps.oeW[3]
    Wran    = ps.oeW[4]
    Wape    = ps.oeW[5]
    Wp      = ps.Wp
    f       = ps.tMax / m
	mu 		= ps.μ
    rpermin = ps.rpmin
    k_petro = ps.k

    # Begin generated code
    t2 = cos(ape);
    t3 = cos(inc);
    t4 = sin(ape);
    t5 = sin(inc);
    t10 = b_petro+1.0;
    t11 = e+1.0;
    t12 = e*e;
    t13 = sma*sma;
    t17 = -ape_t;
    t18 = e-1.0;
    t19 = 1.0/e;
    t26 = -e_t;
    t27 = 1.0/(f*f);
    t28 = -inc_t;
    t29 = 1.0/m_petro;
    t30 = -ran_t;
    t31 = 1.0/r_petro;
    t32 = 1.0/rpermin;
    t33 = -sma_t;
    t34 = 1.0/sma;
    t37 = 1.0/sma_t;
    t6 = fabs(t2);
    t7 = fabs(t3);
    t8 = fabs(t4);
    t9 = (t4/fabs(t4));
    t14 = t2*t2;
    t15 = t4*t4;
    t16 = t5*t5;
    t20 = 1.0/t12;
    t21 = t19*t19*t19;
    t23 = t19*t19*t19*t19*t19;
    t25 = t19*t19*t19*t19*t19*t19*t19;
    t35 = 1.0/t13;
    t36 = t34*t34*t34;
    t38 = t10*t10;
    t39 = 1.0/t5;
    t42 = 1.0/t11;
    t43 = t12-1.0;

    #t44 = ape+t17;
	if (fmod(ape, 2.0*pi) > fmod(ape_t, 2.0*pi))
		t44 = fmod(ape - ape_t, 2.0*pi) >  1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : 1.0e-6;
	else
		t44 = fmod(ape - ape_t, 2.0*pi) < -1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : -1.0e-6;
    end

    t45 = e+t26;
    t46 = inc+t28;

    #t47 = ran+t30;
	if (fmod(ran, 2.0*pi) > fmod(ran_t, 2.0*pi))
		t47 = fmod(ran - ran_t, 2.0*pi) >  1.0e-6 ? fmod(ran - ran_t, 2.0*pi) :  1.0e-6;
	else
		t47 = fmod(ran - ran_t, 2.0*pi) < -1.0e-6 ? fmod(ran - ran_t, 2.0*pi) : -1.0e-6;
    end

    t48 = sma+t33;
    t69 = sma*t18*t32;
    t22 = t20*t20;
    t24 = t20*t20*t20;
    t40 = e*t6;
    t41 = e*t8;
    t49 = cos(t47);
    t50 = -t20;
    t51 = e*t2*t9;
    t52 = cos(t44);
    t55 = t48*t48;
    t57 = t12*t14;
    t58 = t12*t15;
    t60 = t43*t43;
    t61 = t45*t45;
    t62 = t46*t46;
    t63 = sma*t43;
    t65 = 1.0/t43;
    t72 = t23*t43;
    t76 = t69+1.0;
    t77 = (t21*t43)/2.0;
    t82 = t29*t37*t48;
    t53 = -t40;
    t54 = -t41;
    t56 = acos(t52);
    t59 = acos(t49);
    t64 = mu*t63;
    t66 = 1.0/t60;
    t67 = -t57;
    t68 = -t58;
    t78 = t22*t43*(3.0/2.0);
    t79 = k_petro*t76;
    t80 = -t77;
    t86 = (t24*t60)/4.0;
    t87 = t25*t60*(3.0/2.0);
    t92 = pow(t82,n_petro);
    t112 = (We*mu*t27*t34*t61*t65)/4.0;
    t70 = t56*t56;
    t71 = t59*t59;
    t73 = -t64;
    t74 = t67+1.0;
    t75 = t68+1.0;
    t81 = -t78;
    t83 = exp(t79);
    t93 = -t87;
    t96 = t92+1.0;
    t97 = t86+1.0/2.7E+1;
    t84 = sqrt(t74);
    t85 = sqrt(t75);
    t88 = 1.0/sqrt(t73);
    t94 = Wp*t83;
    t99 = pow(t96,t31);
    t102 = sqrt(t97);
    t111 = t72+t93;
    t89 = t88*t88*t88;
    t90 = 1.0/t84;
    t91 = 1.0/t85;
    t95 = t94+1.0;
    t100 = t54+t84;
    t101 = t53+t85;
    t104 = 1.0/t102;
    t105 = pow(t41-t84,2.0);
    t106 = pow(t40-t85,2.0);
    t108 = -1.0/(t41-t84);
    t114 = t77+t102;
    t115 = t80+t102;
    t122 = (Wsma*mu*t18*t27*t36*t42*t55*t99)/4.0;
    t126 = (b_petro*f*t7*t39*t63*t88)/(t41-t84);
    t98 = e*t14*t90;
    t103 = t2*t4*t12*t90;
    t109 = 1.0/t105;
    t116 = pow(t114,1.0/3.0);
    t117 = pow(t115,1.0/3.0);
    t121 = Winc*mu*t27*t34*t62*t65*t106;
    t123 = (t104*t111)/2.0;
    t124 = Wran*mu*t16*t27*t34*t65*t71*t105;
    t125 = b_petro*f*t7*t39*t63*t88*t108;
    t107 = t8+t98;
    t110 = -t103;
    t118 = 1.0/(t116*t116);
    t119 = 1.0/(t117*t117);
    t120 = -t117;
    t127 = t20+t81+t123;
    t128 = t50+t78+t123;
    t113 = t51+t110;
    t129 = t19+t116+t120;
    t140 = (t118*t127)/3.0;
    t142 = (t119*t128)/3.0;
    t130 = e*t129;
    t131 = t129*t129;
    t141 = -t140;
    t132 = t131-1.0;
    t133 = t130-1.0;
    t135 = t13*t60*t131;
    t145 = t20+t141+t142;
    t134 = 1.0/t133;
    t136 = t63*t134;
    t137 = -t136;
    t138 = t63+t137;
    t139 = t138*t138;
    t143 = t132*t139;
    t144 = -t143;
    t146 = t135+t144;
    t147 = sqrt(t146);
    t148 = 1.0/t147;
    t149 = f*t19*t88*t147;
    t150 = t126+t149;
    t151 = 1.0/(t150*t150);
    t152 = 1.0/(t150*t150*t150);
    t153 = Wape*t38*t70*t151;
    t154 = -t153;
    t155 = t112+t121+t122+t124+t154;

    dQdsma  = t95*(Wape*t38*t70*t152*(f*t19*t88*t148*(t132*t138*(t43-t43*t134)*2.0-sma*t60*t131*2.0)*(-1.0/2.0)+(f*mu*t19*t43*t89*t147)/2.0+(b_petro*f*t7*t39*t43*t88)/(t41-t84)+(b_petro*f*mu*sma*t7*t39*t60*t89)/(t41*2.0-t84*2.0))*-2.0+(We*mu*t27*t35*t61*t65)/4.0+Winc*mu*t27*t35*t62*t65*t106+Wsma*mu*t18*t27*(t35*t35)*t42*t55*t99*(3.0/4.0)-(Wsma*mu*t18*t27*t36*t42*t99*(sma*2.0-sma_t*2.0))/4.0+Wran*mu*t16*t27*t35*t65*t71*t105-(Wsma*mu*n_petro*t18*t27*t29*t31*t36*t37*t42*t55*pow(t82,n_petro-1.0)*pow(t96,t31-1.0))/4.0)-k_petro*t18*t32*t94*t155;
    dQde    = -t95*(Wape*t38*t70*t152*(f*t50*t88*t147+(f*t19*t88*t148*(t129*t139*t145*2.0-t132*t138*(e*sma*2.0+t134*t136*(t129-e*t145)-e*sma*t134*2.0)*2.0+t13*t43*t129*t130*4.0-t13*t60*t129*t145*2.0))/2.0+f*mu*sma*t89*t147+(b_petro*e*f*sma*t7*t39*t88*2.0)/(t41-t84)-b_petro*f*t7*t39*t63*t88*t107*t109+(b_petro*e*f*mu*t7*t13*t39*t43*t89)/(t41-t84))*2.0+(We*mu*t27*t34*t65*(e*2.0-e_t*2.0))/4.0-(We*e*mu*t27*t34*t61*t66)/2.0+(Wsma*mu*t27*t36*t42*t55*t99)/4.0+Winc*mu*t27*t34*t62*t65*(t6+e*t15*t91)*(t40-t85)*2.0-(Wsma*mu*t18*t27*t36*(t42*t42)*t55*t99)/4.0-Winc*e*mu*t27*t34*t62*t66*t106*2.0-Wran*e*mu*t16*t27*t34*t66*t71*t105*2.0+Wran*mu*t16*t27*t34*t65*t71*t107*(t41-t84)*2.0)-k_petro*sma*t32*t94*t155;
    dQdinc  = -t95*(Wape*t38*t70*t152*((b_petro*f*t63*t88*((t3/fabs(t3))))/(t41-t84)+(b_petro*f*t3*t7*t63*t88)/(t16*(t41-t84)))*-2.0+Winc*mu*t27*t34*t65*t106*(inc*2.0-inc_t*2.0)+Wran*mu*t3*t5*t27*t34*t65*t71*t105*2.0);
    dQdape  = t95*(Wape*t38*t56*t151*sin(t44)*1.0/sqrt(-t52*t52+1.0)*2.0+Winc*mu*t27*t34*t62*t65*(t40-t85)*(e*t4*((t2/fabs(t2)))-t2*t4*t12*t91)*2.0-Wran*mu*t16*t27*t34*t65*t71*t113*(t41-t84)*2.0+Wape*b_petro*f*t7*t38*t39*t63*t70*t88*t109*t113*t152*2.0);
    dQdran  = Wran*mu*t16*t27*t34*t59*t65*t95*t105*sin(t47)*1.0/sqrt(-t49*t49+1.0)*-2.0;

    return SVector(dQdsma, dQde, dQdinc, dQdran, dQdape, 0.0)
end

