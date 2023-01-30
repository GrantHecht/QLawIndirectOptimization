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
	#t109 = ape - ape_t;
	if (fmod(ape, 2.0*pi) > fmod(ape_t, 2.0*pi))
		t109 = fmod(ape - ape_t, 2.0*pi) >  1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : 1.0e-6;
	else
		t109 = fmod(ape - ape_t, 2.0*pi) < -1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : -1.0e-6;
    end

	t110 = cos(t109);
	t2 = acos(t110);
	t3 = e*e;
	t4 = t3 - 1.0;
	t5 = 1.0 / (e*e*e);
	t6 = 1.0 / (e*e*e*e*e*e);
	t7 = t4*t4;
	t8 = t6*t7*(1.0 / 4.0);
	t9 = t8 + 1.0 / 2.7e1;
	t10 = sqrt(t9);
	t11 = 1.0 / e;
	t13 = t4*t5*(1.0 / 2.0);
	t14 = t10 + t13;
	t15 = pow(t14, 1.0 / 3.0);
	t12 = t11 + t15 - pow(t10 - t4*t5*(1.0 / 2.0), 1.0 / 3.0);
	t17 = t10 - t13;
	t18 = pow(t17, 1.0 / 3.0);
	t19 = t11 + t15 - t18;
	t24 = sma*t4;
	t25 = e*t19;
	t26 = t25 - 1.0;
	t27 = 1.0 / t26;
	t28 = sma*t4*t27;
	t16 = t24 - t28;
	t20 = cos(ape);
	t32 = mu*sma*t4;
	t21 = 1.0 / sqrt(-t32);
	t22 = b_petro + 1.0;
	t23 = t19*t19;
	t29 = t16*t16;
	t30 = sma*sma;
	t31 = t7*t23*t30;
	t33 = t23 - 1.0;
	t36 = t29*t33;
	t34 = t31 - t36;
	t35 = sqrt(t34);
	t37 = 1.0 / (e*e);
	t38 = 1.0 / sqrt(t9);
	t39 = 1.0 / (e*e*e*e*e);
	t40 = t4*t39;
	t41 = 1.0 / (e*e*e*e*e*e*e);
	t46 = t7*t41*(3.0 / 2.0);
	t42 = t40 - t46;
	t43 = t38*t42*(1.0 / 2.0);
	t44 = 1.0 / (e*e*e*e);

	#t45 = 1.0 / pow(t14, 2.0 / 3.0);
	t45 = t14 < 1.0e-10 ? 1.0e20 : 1.0 / pow(t14, 2.0 / 3.0);

	t47 = t4*t44*(3.0 / 2.0);
	t48 = 1.0 / pow(t17, 2.0 / 3.0);
	t49 = -t37 + t43 + t47;
	t50 = t48*t49*(1.0 / 3.0);
	t51 = t37 + t43 - t47;
	t53 = t45*t51*(1.0 / 3.0);
	t52 = t37 + t50 - t53;
	t54 = cos(inc);
	t55 = fabs(t54);
	t56 = sin(inc);
	t57 = 1.0 / t56;
	t58 = t20*t20;
	t65 = t3*t58;
	t59 = -t65 + 1.0;
	t60 = sqrt(t59);
	t61 = sin(ape);
	t62 = fabs(t61);
	t66 = e*t62;
	t63 = t60 - t66;
	t64 = 1.0 / t63;
	t67 = 1.0 / pow(-t32, 3.0 / 2.0);
	t68 = 1.0 / (f*f);
	t69 = sma - sma_t;
	t70 = 1.0 / sma;
	t71 = e - e_t;
	t72 = 1.0 / t4;
	t73 = inc - inc_t;
	t74 = fabs(t20);
	t75 = t61*t61;
	t78 = t3*t75;
	t76 = -t78 + 1.0;
	t77 = 1.0 / (t4*t4);
	t79 = sqrt(t76);
	t114 = e*t74;
	t80 = t79 - t114;
	t81 = t73*t73;
	t82 = 1.0 / (sma*sma*sma);
	t83 = t69*t69;
	t84 = e - 1.0;
	t85 = e + 1.0;
	t86 = fabs(t69);
	t87 = fabs(m_petro);
	t88 = 1.0 / t87;
	t89 = fabs(sma_t);
	t90 = 1.0 / t89;
	t91 = t86*t88*t90;
	t92 = pow(t91, n_petro);
	t93 = t92 + 1.0;
	t94 = 1.0 / r_petro;
	t95 = pow(t93, t94);

	#double t100 = ran-ran_t;
	if (fmod(ran, 2.0*pi) > fmod(ran_t, 2.0*pi))
		t100 = fmod(ran - ran_t, 2.0*pi) >  1.0e-6 ? fmod(ran - ran_t, 2.0*pi) :  1.0e-6;
	else
		t100 = fmod(ran - ran_t, 2.0*pi) < -1.0e-6 ? fmod(ran - ran_t, 2.0*pi) : -1.0e-6;
    end

	t101 = cos(t100);
	t96 = acos(t101);
	t97 = 1.0 / sqrt(t59);
	t98 = e*t58*t97;
	t99 = t62 + t98;
	t102 = t96*t96;
	t103 = t56*t56;
	t104 = 1.0 / rpermin;
	t105 = sma*t84*t104;
	t106 = t105 + 1.0;
	t107 = k_petro*t106;
	t108 = exp(t107);
	t111 = t2*t2;
	t112 = t22*t22;
	t113 = t71*t71;
	t115 = t80*t80;
	t116 = t63*t63;
	t117 = 1.0 / t85;
	t118 = sin(tru);
	t119 = Wp*t108;
	t120 = t119 + 1.0;
	t121 = f*t11*t21*t35;
	t125 = b_petro*f*sma*t4*t21*t55*t57*t64;
	t122 = t121 - t125;
	t123 = 1.0 / sqrt(t34);
	t124 = 1.0 / (sma*sma);
	t126 = 1.0 / (t122*t122);
	t127 = We*mu*t68*t70*t72*t113*(1.0 / 4.0);
	t128 = Winc*mu*t68*t70*t72*t81*t115;
	t129 = Wran*mu*t68*t70*t72*t102*t103*t116;
	t130 = Wsma*mu*t68*t82*t83*t84*t95*t117*(1.0 / 4.0);
	t163 = Wape*t111*t112*t126;
	t131 = t127 + t128 + t129 + t130 - t163;
	t132 = 1.0 / sqrt(t76);
	t133 = 1.0 / (t63*t63);
	t134 = 1.0 / (t122*t122*t122);

	#double t135 = (t61 / fabs(t61));
	t135 = t61 > 0.0 ? 1.0 : -1.0;

	t136 = e*t20*t135;
	t202 = t3*t20*t61*t97;
	t137 = t136 - t202;
	t138 = f*t21*t35*t37;
	t139 = t19*t29*t52*2.0;
	t140 = e*sma*2.0;
	t141 = 1.0 / (t26*t26);
	t142 = t11 + t15 - t18 - e*t52;
	t143 = sma*t4*t141*t142;
	t144 = t140 + t143 - e*sma*t27*2.0;
	t145 = e*t4*t23*t30*4.0;
	t146 = t139 + t145 - t16*t33*t144*2.0 - t7*t19*t30*t52*2.0;
	t147 = b_petro*e*f*sma*t21*t55*t57*t64*2.0;
	t148 = b_petro*f*sma*t4*t21*t55*t57*t99*t133;
	t149 = b_petro*e*f*mu*t4*t30*t55*t57*t64*t67;
	t150 = t138 + t147 + t148 + t149 - f*mu*sma*t35*t67 - f*t11*t21*t123*t146*(1.0 / 2.0);
	t151 = e*2.0;
	t152 = e_t*2.0;
	t153 = t151 - t152;
	t154 = We*e*mu*t68*t70*t77*t113*(1.0 / 2.0);
	t155 = e*t75*t132;
	t156 = t74 + t155;
	t157 = Winc*mu*t68*t70*t72*t80*t81*t156*2.0;
	t158 = Winc*e*mu*t68*t70*t77*t81*t115*2.0;
	t159 = 1.0 / (t85*t85);
	t160 = Wsma*mu*t68*t82*t83*t84*t95*t159*(1.0 / 4.0);
	t161 = Wran*mu*t63*t68*t70*t72*t99*t102*t103*2.0;
	t162 = Wran*e*mu*t68*t70*t77*t102*t103*t116*2.0;
	t164 = cos(tru);
	t165 = e*t164;
	t166 = t165 + 1.0;
	t167 = 1.0 / t166;
	t168 = sma*t7*t23*2.0;
	t169 = t4*t27;
	t170 = -t3 + t169 + 1.0;
	t171 = t16*t33*t170*2.0;
	t172 = t168 + t171;
	t173 = f*t11*t21*t123*t172*(1.0 / 2.0);
	t174 = f*mu*t4*t11*t35*t67*(1.0 / 2.0);
	t175 = t173 + t174 - b_petro*f*t4*t21*t55*t57*t64 - b_petro*f*mu*sma*t7*t55*t57*t64*t67*(1.0 / 2.0);
	t176 = We*mu*t68*t72*t113*t124*(1.0 / 4.0);
	t177 = Winc*mu*t68*t72*t81*t115*t124;
	t178 = Wran*mu*t68*t72*t102*t103*t116*t124;
	t179 = sma*2.0;
	t180 = sma_t*2.0;
	t181 = t179 - t180;
	t182 = 1.0 / (sma*sma*sma*sma);
	t183 = Wsma*mu*t68*t83*t84*t95*t117*t182*(3.0 / 4.0);

	#double t184 = (t69 / fabs(t69));
	t184 = t69 > 0.0 ? 1.0 : -1.0;

	t185 = t94 - 1.0;
	t186 = pow(t93, t185);
	t187 = n_petro - 1.0;
	t188 = pow(t91, t187);
	t189 = t176 + t177 + t178 + t183 - Wape*t111*t112*t134*t175*2.0 - Wsma*mu*t68*t82*t84*t95*t117*t181*(1.0 / 4.0) - Wsma*mu*n_petro*t68*t82*t83*t84*t88*t90*t94*t117*t184*t186*t188*(1.0 / 4.0);
	t190 = t120*t189;
	t191 = t190 - Wp*k_petro*t84*t104*t108*t131;
	t192 = sma*t4*t167;
	t193 = t24 + t192;
	t194 = sin(t109);
	t195 = t110*t110;
	t196 = -t195 + 1.0 > 1e-6 ? -t195 + 1.0 : 1e-7;
	t197 = 1.0 / sqrt(t196);
	t198 = Wape*t2*t112*t126*t194*t197*2.0;

	#double t199 = (t20 / fabs(t20));
	t199 = t20 > 0.0 ? 1.0 : -1.0;

	t200 = e*t61*t199;
	t207 = t3*t20*t61*t132;
	t201 = t200 - t207;
	t203 = Wran*mu*t63*t68*t70*t72*t102*t103*t137*2.0;
	t204 = Wape*b_petro*f*sma*t4*t21*t55*t57*t111*t112*t133*t134*t137*2.0;
	t208 = Winc*mu*t68*t70*t72*t80*t81*t201*2.0;
	t205 = t198 + t203 + t204 - t208;
	t206 = ape + tru;
	t209 = sin(t206);

	fr = sma*t4*t21*t118*(t120*(t154 + t157 + t158 + t160 + t161 + t162 + Wape*t111*t112*t150*1.0 / pow(f*t11*t21*sqrt(t31 - t29*(t12*t12 - 1.0)) - b_petro*f*sma*t4*t21*t55*t57*t64, 3.0)*2.0 - We*mu*t68*t70*t72*t153*(1.0 / 4.0) - Wsma*mu*t68*t82*t83*t95*t117*(1.0 / 4.0)) - Wp*k_petro*sma*t104*t108*t131) - e*t21*t30*t118*t191*2.0 - sma*t4*t11*t21*t120*t164*t205;
    fθ = t21*(t120*(t154 + t157 + t158 + t160 + t161 + t162 + Wape*t111*t112*t134*t150*2.0 - We*mu*t68*t70*t72*t153*(1.0 / 4.0) - Wsma*mu*t68*t82*t83*t95*t117*(1.0 / 4.0)) - Wp*k_petro*sma*t104*t108*t131)*(t164*t193 + e*sma*t4*t167) - t21*t30*t166*t191*2.0 + t11*t21*t118*t120*t193*t205;
	fh = -sma*t4*t21*t120*t167*cos(t206)*(Wape*t111*t112*t134*(b_petro*f*sma*t4*t21*t64*((t54 / fabs(t54))) + b_petro*f*sma*t4*t21*t54*t55*1.0 / (t56*t56)*t64)*2.0 + Winc*mu*t68*t70*t72*t115*(inc*2.0 - inc_t*2.0) + Wran*mu*t54*t56*t68*t70*t72*t102*t116*2.0) - sma*t4*t21*t54*t57*t120*t167*t205*t209 - Wran*mu*t21*t56*t68*t96*t116*t120*t167*t209*sin(t100)*1.0 / sqrt(-t101*t101 + 1.0)*2.0;

    α  = atan(fr, fθ)
    β  = atan(fh / sqrt(fr*fr + fθ*fθ))
    return (α, β)
end