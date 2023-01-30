
function  dQn(sma, e, inc, ape, ran, tru, m, alpha, beta, ps::qLawParams)
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
    t3 = cos(beta);
    t4 = cos(inc);
    t5 = sin(ape);
    t6 = cos(tru);
    t7 = sin(inc);
    t8 = sin(tru);
    t11 = ape+tru;
    t15 = b_petro+1.0;
    t16 = e*2.0;
    t17 = e+1.0;
    t18 = e*e;
    t19 = e_t*2.0;
    t20 = sma*2.0;
    t21 = sma_t*2.0;
    t22 = sma*sma;
    t29 = -ape_t;
    t30 = e-1.0;
    t31 = 1.0/e;
    t38 = -e_t;
    t40 = 1.0/(f*f);
    t41 = -inc_t;
    t42 = 1.0/m_petro;
    t43 = n_petro-1.0;
    t44 = -ran_t;
    t45 = 1.0/r_petro;
    t46 = 1.0/rpermin;
    t47 = -sma_t;
    t49 = 1.0/sma;
    t53 = 1.0/sma_t;
    t55 = e*sma*-2.0;
    t9 = fabs(t2);
    t10 = fabs(t4);
    t12 = fabs(t5);

    #t13 = (t2/fabs(t2));
    #t14 = (t5/fabs(t5));
    t13 = t2 > 0.0 ? 1.0 : -1.0
    t14 = t5 > 0.0 ? 1.0 : -1.0

    t23 = t2*t2;
    t24 = t5*t5;
    t25 = sma*t16;
    t26 = t7*t7;
    t27 = e*t6;
    t28 = sin(t11);
    t32 = 1.0/t18;
    t33 = t31*t31*t31;
    t35 = t31*t31*t31*t31*t31;
    t37 = t31*t31*t31*t31*t31*t31*t31;
    t39 = -t19;
    t48 = -t21;
    t50 = 1.0/t22;
    t51 = t49*t49*t49;
    t54 = t15*t15;
    t56 = 1.0/t7;
    t59 = -t18;
    t60 = 1.0/t17;
    t61 = t18-1.0;

    #t63 = ape+t29;
	if (mod(ape, 2.0*pi) > mod(ps.oet[5], 2.0*pi))
		t63 = mod(ape - ps.oet[5], 2.0*pi) >  1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : 1.0e-6;
	else
		t63 = mod(ape - ps.oet[5], 2.0*pi) < -1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : -1.0e-6;
	end
    t63 = acos(cos(t63))

    t64 = e+t38;
    t65 = inc+t41;

    #t66 = ran+t44;
	if (mod(ran, 2.0*pi) > mod(ps.oet[4], 2.0*pi))
		t66 = mod(ran - ps.oet[4], 2.0*pi) >  1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) :  1.0e-6;
	else
		t66 = mod(ran - ps.oet[4], 2.0*pi) < -1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) : -1.0e-6;
	end
    t66 = acos(cos(t66))

    t67 = sma+t47;
    t70 = -t31;
    t72 = t45-1.0;
    t97 = sma*t30*t46;
    t34 = t32*t32;
    t36 = t32*t32*t32;
    t52 = t50*t50;
    t57 = e*t9;
    t58 = e*t12;
    t62 = t60*t60;
    t68 = t27+1.0;
    t69 = cos(t66);
    t71 = -t32;
    t73 = e*t5*t13;
    t74 = e*t2*t14;
    t75 = cos(t63);
    t77 = sin(t63);
    t79 = t20+t48;
    t80 = t67*t67;
    t82 = t18*t23;
    t83 = t18*t24;
    t85 = t61*t61;
    t86 = t16+t39;
    t87 = t64*t64;
    t88 = t65*t65;
    t89 = sma*t61;
    t93 = 1.0/t61;
    t95 = t23*t59;
    t96 = t24*t59;
    t100 = t35*t61;
    t105 = t97+1.0;
    t107 = (t33*t61)/2.0;
    t112 = t42*t53*t67;
    t76 = -t57;
    t78 = -t58;
    t81 = acos(t75);
    t84 = acos(t69);
    t90 = t75*t75;
    t91 = 1.0/t68;
    t92 = mu*t89;
    t94 = 1.0/t85;
    t103 = t95+1.0;
    t104 = t96+1.0;
    t108 = t34*t61*(3.0/2.0);
    t109 = k_petro*t105;
    t110 = -t107;
    t116 = (t36*t85)/4.0;
    t117 = t37*t85*(3.0/2.0);
    t123 = pow(t112,n_petro);
    t129 = pow(t112,t43);
    t151 = (We*mu*t40*t49*t86*t93)/4.0;
    t152 = (We*mu*t40*t49*t87*t93)/4.0;
    t153 = (We*mu*t40*t50*t87*t93)/4.0;
    t98 = t81*t81;
    t99 = t84*t84;
    t101 = -t92;
    t102 = -t90;
    t111 = -t108;
    t113 = exp(t109);

    #t114 = sqrt(t103);
    t114 = t103 < 0.0 ? 0.0 : sqrt(t103)

    #t115 = sqrt(t104);
    t115 = t104 < 0.0 ? 0.0 : sqrt(t104)

    t122 = t89*t91;
    t124 = -t117;
    t128 = t123+1.0;
    t130 = t116+1.0/2.7E+1;
    t154 = -t151;
    t155 = (We*e*mu*t40*t49*t87*t94)/2.0;
    t106 = t102+1.0;
    t118 = 1.0/sqrt(t101);
    t120 = 1.0/t114;
    t121 = 1.0/t115;
    t125 = Wp*t113;
    t133 = pow(t128,t45);
    t134 = t78+t114;
    t135 = t76+t115;
    t136 = sqrt(t130);
    t140 = t89+t122;
    t141 = pow(t58-t114,2.0);
    t142 = pow(t57-t115,2.0);
    t143 = pow(t128,t72);
    t146 = -1.0/(t58-t114);
    t150 = t100+t124;
    t119 = t118*t118*t118;
    t126 = t125+1.0;
    t127 = 1.0/sqrt(t106);
    t131 = e*t24*t121;
    t132 = e*t23*t120;
    t137 = t2*t5*t18*t120;
    t138 = t2*t5*t18*t121;
    t139 = 1.0/t136;
    t147 = 1.0/t141;
    t148 = t2*t5*t59*t120;
    t149 = t2*t5*t59*t121;
    t158 = t107+t136;
    t159 = t110+t136;
    t166 = (Wsma*mu*t40*t51*t60*t80*t133)/4.0;
    t168 = Winc*mu*t40*t49*t88*t93*t142;
    t169 = Winc*mu*t40*t50*t88*t93*t142;
    t170 = Winc*mu*t16*t40*t49*t88*t94*t142;
    t171 = (Wsma*mu*t30*t40*t51*t60*t79*t133)/4.0;
    t173 = (Wsma*mu*t30*t40*t51*t62*t80*t133)/4.0;
    t174 = Wsma*mu*t30*t40*t52*t60*t80*t133*(3.0/4.0);
    t176 = (b_petro*e*f*sma*t10*t56*t118*-2.0)/(t58-t114);
    t178 = Wran*mu*t26*t40*t49*t93*t99*t141;
    t179 = Wran*mu*t26*t40*t50*t93*t99*t141;
    t180 = b_petro*f*t10*t56*t61*t118*t146;
    t181 = b_petro*f*t10*t56*t89*t118*t146;
    t182 = (b_petro*f*t10*t56*t61*t118)/(t58-t114);
    t183 = Wran*mu*t16*t26*t40*t49*t94*t99*t141;
    t184 = (b_petro*f*t10*t56*t89*t118)/(t58-t114);
    t195 = (Wsma*mu*n_petro*t30*t40*t42*t45*t51*t53*t60*t80*t129*t143)/4.0;
    t144 = t12+t132;
    t145 = t9+t131;
    t156 = t73+t149;
    t157 = t74+t148;

    t158 = t158 < 0.0 ? 0.0 : t158
    t160 = pow(t158,1.0/3.0);

    t159 = t159 < 0.0 ? 0.0 : t159
    t161 = pow(t159,1.0/3.0);

    t167 = -t166;
    t172 = t30*t166;
    t175 = -t171;
    t177 = (t139*t150)/2.0;
    t185 = b_petro*e*f*mu*t10*t22*t56*t61*t119*t146;
    t186 = (b_petro*f*mu*sma*t10*t56*t85*t119*(-1.0/2.0))/(t58-t114);
    t187 = (b_petro*f*mu*sma*t10*t56*t85*t119)/(t58*2.0-t114*2.0);
    t196 = -t195;

    #t162 = 1.0/(t160*t160);
    t162 = t160*t160 < 1e-10 ? 1e20 : 1.0 / (t160*t160)

    #t163 = 1.0/(t161*t161);
    t163 = t161*t161 < 1e-10 ? 1e20 : 1.0 / (t161*t161)

    t164 = -t160;
    t165 = -t161;
    t188 = t32+t111+t177;
    t189 = Winc*mu*t40*t49*t88*t93*t145*(t57-t115)*-2.0;
    t190 = t71+t108+t177;
    t191 = Wran*mu*t26*t40*t49*t93*t99*t144*(t58-t114)*-2.0;
    t192 = b_petro*f*t10*t56*t89*t118*t144*t147;
    t193 = Winc*mu*t40*t49*t88*t93*t156*(t57-t115)*-2.0;
    t194 = Winc*mu*t40*t49*t88*t93*t156*(t57-t115)*2.0;
    t197 = Wran*mu*t26*t40*t49*t93*t99*t157*(t58-t114)*-2.0;
    t198 = t31+t160+t165;
    t215 = (t162*t188)/3.0;
    t217 = (t163*t190)/3.0;
    t199 = e*t198;
    t200 = t198*t198;
    t216 = -t215;
    t201 = t200-1.0;
    t202 = t199-1.0;
    t206 = t20*t85*t200;
    t207 = t22*t85*t200;
    t208 = t22*t61*t198*t199*4.0;
    t220 = t32+t216+t217;
    t203 = 1.0/t202;
    t221 = e*t220;
    t231 = sma*t20*t85*t198*t220;
    t232 = t22*t85*t198*t220*-2.0;
    t204 = t203*t203;
    t205 = t25*t203;
    t209 = t61*t203;
    t210 = t89*t203;
    t226 = t70+t161+t164+t221;
    t211 = -t210;
    t212 = t59+t209+1.0;
    t241 = -t89*t204*(t198-t221);
    t213 = t89+t211;
    t244 = t55+t205+t241;
    t214 = t213*t213;
    t225 = t201*t212*t213*2.0;
    t250 = t201*t213*(-t205+e*sma*2.0+t89*t204*(t198-t221))*-2.0;
    t218 = t201*t214;
    t239 = t206+t225;
    t242 = t198*t214*t220*2.0;
    t219 = -t218;
    t256 = t208+t232+t242+t250;
    t222 = t207+t219;

    #t223 = sqrt(t222);
    t223 = t222 < 0.0 ? 0.0 : sqrt(t222)

    t224 = 1.0/t223;
    t227 = f*mu*sma*t119*t223;
    t228 = f*t31*t118*t223;
    t229 = f*t32*t118*t223;
    t233 = (f*mu*t31*t61*t119*t223)/2.0;
    t230 = -t227;
    t234 = t184+t228;
    t251 = (f*t31*t118*t224*t239)/2.0;
    t259 = (f*t31*t118*t224*t256)/2.0;
    t235 = 1.0/(t234*t234);
    t236 = 1.0/(t234*t234*t234);
    t253 = t182+t187+t233+t251;
    t260 = -t259;
    t237 = Wape*t54*t98*t235;
    t240 = Wape*t54*t77*t81*t127*t235*2.0;
    t243 = Wape*b_petro*f*t10*t20*t54*t56*t61*t98*t118*t147*t157*t236;
    t254 = Wape*t54*t98*t236*t253*2.0;
    t262 = t176+t185+t192+t229+t230+t260;
    t263 = Wape*t54*t98*t236*(t185+t192+t229+t230+t260+b_petro*f*t10*t25*t56*t118*t146)*2.0;
    t238 = -t237;
    t252 = t194+t197+t240+t243;
    t255 = -t254;
    t264 = t154+t155+t167+t170+t173+t183+t189+t191+t263;
    t265 = t126*(t154+t155+t167+t173+t189+t191+t263+Winc*e*mu*t40*t49*t88*t94*t142*2.0+Wran*e*mu*t26*t40*t49*t94*t99*t141*2.0);
    t245 = t152+t168+t172+t178+t238;
    t257 = t153+t169+t174+t175+t179+t196+t255;
    t246 = k_petro*sma*t46*t125*t245;
    t248 = k_petro*t30*t46*t125*t245;
    t258 = t126*t257;
    t247 = -t246;
    t249 = -t248;
    val = sin(beta)*(t118*t122*t126*cos(t11)*(Wape*t54*t98*t236*((b_petro*f*t89*t118*((t4/fabs(t4))))/(t58-t114)+(b_petro*f*t4*t10*t89*t118)/(t26*(t58-t114)))*-2.0+Winc*mu*t40*t49*t93*t142*(inc*2.0-inc_t*2.0)+Wran*mu*t4*t7*t40*t49*t93*t99*t141*2.0)+t4*t28*t56*t118*t122*t126*t252+Wran*mu*t7*t28*t40*t84*t91*t118*t126*t141*sin(t66)*1.0/sqrt(-t69*t69+1.0)*2.0)-t3*cos(alpha)*(-t118*(t246-t265)*(e*t122+t6*t140)+sma*t20*t68*t118*(t248-t258)+t8*t31*t118*t126*t140*t252)+t3*sin(alpha)*(t8*t89*t118*(t246-t265)-t8*t16*t22*t118*(t248-t258)+t6*t31*t89*t118*t126*t252);

    return val
end