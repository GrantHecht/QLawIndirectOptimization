
function qLawThrustAngles(sma, e, inc, ape, ran, tru, m, ps::qLawParams)
    # Set some constants
    b_petro = 0.01
    m_petro = 3
    n_petro = 4
    r_petro = 2

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

    # Start of generated code
    t2 = cos(ape);
    t3 = cos(inc);
    t4 = sin(ape);
    t5 = cos(tru);
    t6 = sin(inc);
    t7 = sin(tru);
    t10 = ape+tru;
    t14 = b_petro+1.0;
    t15 = e*2.0;
    t16 = e+1.0;
    t17 = e*e;
    t18 = e_t*2.0;
    t19 = sma*2.0;
    t20 = sma_t*2.0;
    t21 = sma*sma;
    t28 = -ape_t;
    t29 = e-1.0;
    t30 = 1.0/e;
    t37 = -e_t;
    t39 = 1.0/(f*f);
    t40 = -inc_t;
    t41 = 1.0/m_petro;
    t42 = n_petro-1.0;
    t43 = -ran_t;
    t44 = 1.0/r_petro;
    t45 = 1.0/rpermin;
    t46 = -sma_t;
    t48 = 1.0/sma;
    t52 = 1.0/sma_t;
    t54 = e*sma*-2.0;
    t8 = fabs(t2);
    t9 = fabs(t3);
    t11 = fabs(t4);

    #t12 = (t2/fabs(t2));
    #t13 = (t4/fabs(t4));
    t12 = t2 > 0.0 ? 1.0 : -1.0
    t13 = t4 > 0.0 ? 1.0 : -1.0

    t22 = t2*t2;
    t23 = t4*t4;
    t24 = sma*t15;
    t25 = t6*t6;
    t26 = e*t5;
    t27 = sin(t10);
    t31 = 1.0/t17;
    t32 = t30*t30*t30;
    t34 = t30*t30*t30*t30*t30;
    t36 = t30*t30*t30*t30*t30*t30*t30;
    t38 = -t18;
    t47 = -t20;
    t49 = 1.0/t21;
    t50 = t48*t48*t48;
    t53 = t14*t14;
    t55 = 1.0/t6;
    t58 = -t17;
    t59 = 1.0/t16;
    t60 = t17-1.0;

    #t62 = ape+t28;
	if (mod(ape, 2.0*pi) > mod(ps.oet[5], 2.0*pi))
	    t62 = mod(ape - ps.oet[5], 2.0*pi) >  1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : 1.0e-6;
	else
		t62 = mod(ape - ps.oet[5], 2.0*pi) < -1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : -1.0e-6;
	end

    t63 = e+t37;
    t64 = inc+t40;

    #t65 = ran+t43;
	if (mod(ran, 2.0*pi) > mod(ps.oet[4], 2.0*pi))
	    t65 = mod(ran - ps.oet[4], 2.0*pi) >  1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) :  1.0e-6;
	else
		t65 = mod(ran - ps.oet[4], 2.0*pi) < -1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) : -1.0e-6;
	end

    t66 = sma+t46;
    t69 = -t30;
    t71 = t44-1.0;
    t96 = sma*t29*t45;
    t33 = t31*t31;
    t35 = t31*t31*t31;
    t51 = t49*t49;
    t56 = e*t8;
    t57 = e*t11;
    t61 = t59*t59;
    t67 = t26+1.0;
    t68 = cos(t65);
    t70 = -t31;
    t72 = e*t4*t12;
    t73 = e*t2*t13;
    t74 = cos(t62);
    t76 = sin(t62);
    t78 = t19+t47;
    t79 = t66*t66;
    t81 = t17*t22;
    t82 = t17*t23;
    t84 = t60*t60;
    t85 = t15+t38;
    t86 = t63*t63;
    t87 = t64*t64;
    t88 = sma*t60;
    t92 = 1.0/t60;
    t94 = t22*t58;
    t95 = t23*t58;
    t99 = t34*t60;
    t104 = t96+1.0;
    t106 = (t32*t60)/2.0;
    t111 = t41*t52*t66;
    t75 = -t56;
    t77 = -t57;
    t80 = acos(t74);
    t83 = acos(t68);
    t89 = t74*t74;
    t90 = 1.0/t67;
    t91 = mu*t88;
    t93 = 1.0/t84;
    t102 = t94+1.0;
    t103 = t95+1.0;
    t107 = t33*t60*(3.0/2.0);
    t108 = k_petro*t104;
    t109 = -t106;
    t115 = (t35*t84)/4.0;
    t116 = t36*t84*(3.0/2.0);
    t122 = pow(t111,n_petro);
    t128 = pow(t111,t42);
    t150 = (We*mu*t39*t48*t85*t92)/4.0;
    t151 = (We*mu*t39*t48*t86*t92)/4.0;
    t152 = (We*mu*t39*t49*t86*t92)/4.0;
    t97 = t80*t80;
    t98 = t83*t83;
    t100 = -t91;
    t101 = -t89;
    t110 = -t107;
    t112 = exp(t108);

    #t113 = sqrt(t102);
    t113 = t102 < 0.0 ? 0.0 : sqrt(t102);

    #t114 = sqrt(t103);
    t114 = t103 < 0.0 ? 0.0 : sqrt(t103)

    t121 = t88*t90;
    t123 = -t116;
    t127 = t122+1.0;
    t129 = t115+1.0/2.7e1;
    t153 = -t150;
    t154 = (We*e*mu*t39*t48*t86*t93)/2.0;
    t105 = t101+1.0;
    t117 = 1.0/sqrt(t100);
    t119 = 1.0/t113;
    t120 = 1.0/t114;
    t124 = Wp*t112;
    t132 = pow(t127,t44);
    t133 = t77+t113;
    t134 = t75+t114;
    t135 = sqrt(t129);
    t139 = t88+t121;
    t140 = pow(t57-t113,2.0);
    t141 = pow(t56-t114,2.0);
    t142 = pow(t127,t71);
    t145 = -1.0/(t57-t113);
    t149 = t99+t123;
    t118 = t117*t117*t117;
    t125 = t124+1.0;
    t126 = 1.0/sqrt(t105);
    t130 = e*t23*t120;
    t131 = e*t22*t119;
    t136 = t2*t4*t17*t119;
    t137 = t2*t4*t17*t120;
    t138 = 1.0/t135;
    t146 = 1.0/t140;
    t147 = t2*t4*t58*t119;
    t148 = t2*t4*t58*t120;
    t157 = t106+t135;
    t158 = t109+t135;
    t165 = (Wsma*mu*t39*t50*t59*t79*t132)/4.0;
    t167 = Winc*mu*t39*t48*t87*t92*t141;
    t168 = Winc*mu*t39*t49*t87*t92*t141;
    t169 = Winc*mu*t15*t39*t48*t87*t93*t141;
    t170 = (Wsma*mu*t29*t39*t50*t59*t78*t132)/4.0;
    t172 = (Wsma*mu*t29*t39*t50*t61*t79*t132)/4.0;
    t173 = Wsma*mu*t29*t39*t51*t59*t79*t132*(3.0/4.0);
    t175 = (b_petro*e*f*sma*t9*t55*t117*-2.0)/(t57-t113);
    t177 = Wran*mu*t25*t39*t48*t92*t98*t140;
    t178 = Wran*mu*t25*t39*t49*t92*t98*t140;
    t179 = b_petro*f*t9*t55*t60*t117*t145;
    t180 = b_petro*f*t9*t55*t88*t117*t145;
    t181 = (b_petro*f*t9*t55*t60*t117)/(t57-t113);
    t182 = Wran*mu*t15*t25*t39*t48*t93*t98*t140;
    t183 = (b_petro*f*t9*t55*t88*t117)/(t57-t113);
    t194 = (Wsma*mu*n_petro*t29*t39*t41*t44*t50*t52*t59*t79*t128*t142)/4.0;
    t143 = t11+t131;
    t144 = t8+t130;
    t155 = t72+t148;
    t156 = t73+t147;

    t157 = t157 < 0.0 ? 0.0 : t157
    t159 = pow(t157,1.0/3.0);

    t158 = t158 < 0.0 ? 0.0 : t158
    t160 = pow(t158,1.0/3.0);

    t166 = -t165;
    t171 = t29*t165;
    t174 = -t170;
    t176 = (t138*t149)/2.0;
    t184 = b_petro*e*f*mu*t9*t21*t55*t60*t118*t145;
    t185 = (b_petro*f*mu*sma*t9*t55*t84*t118*(-1.0/2.0))/(t57-t113);
    t186 = (b_petro*f*mu*sma*t9*t55*t84*t118)/(t57*2.0-t113*2.0);
    t195 = -t194;

    #t161 = 1.0/(t159*t159);
    t161 = t159*t159 < 1e-10 ? 1e20 : 1.0 / (t159*t159)

    #t162 = 1.0/(t160*t160);
    t162 = t160*t160 < 1e-10 ? 1e20 : 1.0 / (t160*t160)

    t163 = -t159;
    t164 = -t160;
    t187 = t31+t110+t176;
    t188 = Winc*mu*t39*t48*t87*t92*t144*(t56-t114)*-2.0;
    t189 = t70+t107+t176;
    t190 = Wran*mu*t25*t39*t48*t92*t98*t143*(t57-t113)*-2.0;
    t191 = b_petro*f*t9*t55*t88*t117*t143*t146;
    t192 = Winc*mu*t39*t48*t87*t92*t155*(t56-t114)*-2.0;
    t193 = Winc*mu*t39*t48*t87*t92*t155*(t56-t114)*2.0;
    t196 = Wran*mu*t25*t39*t48*t92*t98*t156*(t57-t113)*-2.0;
    t197 = t30+t159+t164;
    t214 = (t161*t187)/3.0;
    t216 = (t162*t189)/3.0;
    t198 = e*t197;
    t199 = t197*t197;
    t215 = -t214;
    t200 = t199-1.0;
    t201 = t198-1.0;
    t205 = t19*t84*t199;
    t206 = t21*t84*t199;
    t207 = t21*t60*t197*t198*4.0;
    t219 = t31+t215+t216;
    t202 = 1.0/t201;
    t220 = e*t219;
    t230 = sma*t19*t84*t197*t219;
    t231 = t21*t84*t197*t219*-2.0;
    t203 = t202*t202;
    t204 = t24*t202;
    t208 = t60*t202;
    t209 = t88*t202;
    t225 = t69+t160+t163+t220;
    t210 = -t209;
    t211 = t58+t208+1.0;
    t240 = -t88*t203*(t197-t220);
    t212 = t88+t210;
    t243 = t54+t204+t240;
    t213 = t212*t212;
    t224 = t200*t211*t212*2.0;
    t249 = t200*t212*(-t204+e*sma*2.0+t88*t203*(t197-t220))*-2.0;
    t217 = t200*t213;
    t238 = t205+t224;
    t241 = t197*t213*t219*2.0;
    t218 = -t217;
    t255 = t207+t231+t241+t249;
    t221 = t206+t218;

    #t222 = sqrt(t221);
    t222 = t221 < 0.0 ? 0.0 : sqrt(t221)

    t223 = 1.0/t222;
    t226 = f*mu*sma*t118*t222;
    t227 = f*t30*t117*t222;
    t228 = f*t31*t117*t222;
    t232 = (f*mu*t30*t60*t118*t222)/2.0;
    t229 = -t226;
    t233 = t183+t227;
    t250 = (f*t30*t117*t223*t238)/2.0;
    t258 = (f*t30*t117*t223*t255)/2.0;
    t234 = 1.0/(t233*t233);
    t235 = 1.0/(t233*t233*t233);
    t252 = t181+t186+t232+t250;
    t259 = -t258;
    t236 = Wape*t53*t97*t234;
    t239 = Wape*t53*t76*t80*t126*t234*2.0;
    t242 = Wape*b_petro*f*t9*t19*t53*t55*t60*t97*t117*t146*t156*t235;
    t253 = Wape*t53*t97*t235*t252*2.0;
    t261 = t175+t184+t191+t228+t229+t259;
    t262 = Wape*t53*t97*t235*(t184+t191+t228+t229+t259+b_petro*f*t9*t24*t55*t117*t145)*2.0;
    t237 = -t236;
    t251 = t193+t196+t239+t242;
    t254 = -t253;
    t263 = t153+t154+t166+t169+t172+t182+t188+t190+t262;
    t264 = t125*(t153+t154+t166+t172+t188+t190+t262+Winc*e*mu*t39*t48*t87*t93*t141*2.0+Wran*e*mu*t25*t39*t48*t93*t98*t140*2.0);
    t244 = t151+t167+t171+t177+t237;
    t256 = t152+t168+t173+t174+t178+t195+t254;
    t245 = k_petro*sma*t45*t124*t244;
    t247 = k_petro*t29*t45*t124*t244;
    t257 = t125*t256;
    t246 = -t245;
    t248 = -t247;
    t600 = t3 > 0.0 ? 1.0 : -1.0
    D2 = t7*t88*t117*(t245-t264)-t7*t15*t21*t117*(t247-t257)+t5*t30*t88*t117*t125*t251;
    D1 = t117*(t245-t264)*(e*t121+t5*t139)-sma*t19*t67*t117*(t247-t257)+t7*t69*t117*t125*t139*t251;
    D3 = t117*t121*t125*cos(t10)*(Wape*t53*t97*t235*((b_petro*f*t88*t117*(t600))/(t57-t113)+(b_petro*f*t3*t9*t88*t117)/(t25*(t57-t113)))*-2.0+Winc*mu*t39*t48*t92*t141*(inc*2.0-inc_t*2.0)+Wran*mu*t3*t6*t39*t48*t92*t98*t140*2.0)+t3*t27*t55*t117*t121*t125*t251+Wran*mu*t6*t27*t39*t83*t90*t117*t125*t140*sin(t65)*1.0/sqrt(-t68*t68+1.0)*2.0;

    α  = atan(-D2, -D1)
    β  = atan(-D3 / sqrt(D1*D1 + D2*D2))
    return (α, β)
end