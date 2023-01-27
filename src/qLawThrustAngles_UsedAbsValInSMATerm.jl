
function qLawThrustAngles(sma, e, inc, ape, ran, tru, ps::qLawParams)
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
    f       = ps.tMax
	mu 		= ps.μ
    rpermin = ps.rpmin
    k_petro = ps.k

    # Start of generated code
    t2 = cos(ape);
    t3 = cos(inc);
    t4 = fabs(m_petro);
    t5 = fabs(sma_t);
    t6 = sin(ape);
    t7 = cos(tru);
    t8 = sin(inc);
    t9 = sin(tru);
    t12 = ape+tru;
    t16 = b_petro+1.0;
    t17 = e*2.0;
    t18 = e+1.0;
    t19 = e*e;
    t20 = e_t*2.0;
    t21 = sma*2.0;
    t22 = sma_t*2.0;
    t23 = sma*sma;
    t30 = -ape_t;
    t31 = e-1.0;
    t32 = 1.0/e;
    t39 = -e_t;
    t41 = 1.0/(f*f);
    t42 = -inc_t;
    t43 = n_petro-1.0;
    t44 = -ran_t;
    t45 = 1.0/r_petro;
    t46 = 1.0/rpermin;
    t47 = -sma_t;
    t49 = 1.0/sma;
    t56 = e*sma*-2.0;
    t10 = fabs(t2);
    t11 = fabs(t3);
    t13 = fabs(t6);

    #t14 = (t2/fabs(t2));
    t14 = t2 > 0.0 ? 1.0 : -1.0

    #t15 = (t6/fabs(t6));
    t15 = t6 > 0.0 ? 1.0 : -1.0

    t24 = t2*t2;
    t25 = t6*t6;
    t26 = sma*t17;
    t27 = t8*t8;
    t28 = e*t7;
    t29 = sin(t12);
    t33 = 1.0/t19;
    t34 = t32*t32*t32;
    t36 = t32*t32*t32*t32*t32;
    t38 = t32*t32*t32*t32*t32*t32*t32;
    t40 = -t20;
    t48 = -t22;
    t50 = 1.0/t23;
    t51 = t49*t49*t49;
    t53 = t16*t16;
    t54 = 1.0/t4;
    t55 = 1.0/t5;
    t57 = 1.0/t8;
    t60 = -t19;
    t61 = 1.0/t18;
    t62 = t19-1.0;

    #t64 = ape + t30
    #t64 = ape - ape_t;
	if (mod(ape, 2.0*pi) > mod(ps.oet[5], 2.0*pi))
		t64 = mod(ape - ps.oet[5], 2.0*pi) >  1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : 1.0e-6;
	else
		t64 = mod(ape - ps.oet[5], 2.0*pi) < -1.0e-6 ? mod(ape - ps.oet[5], 2.0*pi) : -1.0e-6;
	end

    t65 = e+t39;
    t66 = inc+t42;

    #t67 = ran+t44;
    #t67 = ran - ran_t
	if (mod(ran, 2.0*pi) > mod(ps.oet[4], 2.0*pi))
		t67 = mod(ran - ps.oet[4], 2.0*pi) >  1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) :  1.0e-6;
	else
		t67 = mod(ran - ps.oet[4], 2.0*pi) < -1.0e-6 ? mod(ran - ps.oet[4], 2.0*pi) : -1.0e-6;
	end

    t68 = sma+t47;
    t73 = -t32;
    t75 = t45-1.0;
    t100 = sma*t31*t46;
    t35 = t33*t33;
    t37 = t33*t33*t33;
    t52 = t50*t50;
    t58 = e*t10;
    t59 = e*t13;
    t63 = t61*t61;
    t69 = t28+1.0;
    t70 = cos(t67);
    t71 = fabs(t68);

    #t72 = (t68/fabs(t68));
    t72 = t68 > 0.0 ? 1.0 : -1.0

    t74 = -t33;
    t76 = e*t6*t14;
    t77 = e*t2*t15;
    t78 = cos(t64);
    t80 = sin(t64);
    t82 = t21+t48;
    t83 = t68*t68;
    t85 = t19*t24;
    t86 = t19*t25;
    t88 = t62*t62;
    t89 = t17+t40;
    t90 = t65*t65;
    t91 = t66*t66;
    t92 = sma*t62;
    t96 = 1.0/t62;
    t98 = t24*t60;
    t99 = t25*t60;
    t103 = t36*t62;
    t108 = t100+1.0;
    t110 = (t34*t62)/2.0;
    t79 = -t58;
    t81 = -t59;
    t84 = acos(t78);
    t87 = acos(t70);
    t93 = t78*t78;
    t94 = 1.0/t69;
    t95 = mu*t92;
    t97 = 1.0/t88;
    t106 = t98+1.0;
    t107 = t99+1.0;
    t111 = t35*t62*(3.0/2.0);
    t112 = k_petro*t108;
    t113 = -t110;
    t118 = (t37*t88)/4.0;
    t119 = t38*t88*(3.0/2.0);
    t127 = t54*t55*t71;
    t154 = (We*mu*t41*t49*t89*t96)/4.0;
    t155 = (We*mu*t41*t49*t90*t96)/4.0;
    t156 = (We*mu*t41*t50*t90*t96)/4.0;
    t101 = t84*t84;
    t102 = t87*t87;
    t104 = -t95;
    t105 = -t93;
    t114 = -t111;
    t115 = exp(t112);
    t116 = sqrt(t106);
    t117 = sqrt(t107);
    t124 = t92*t94;
    t125 = -t119;
    t130 = pow(t127,n_petro);
    t131 = t118+1.0/2.7e1;
    t135 = pow(t127,t43);
    t157 = -t154;
    t158 = (We*e*mu*t41*t49*t90*t97)/2.0;
    t109 = t105+1.0;
    t120 = 1.0/sqrt(t104);
    t122 = 1.0/t116;
    t123 = 1.0/t117;
    t126 = Wp*t115;
    t132 = t130+1.0;
    t136 = t81+t116;
    t137 = t79+t117;
    t138 = sqrt(t131);
    t142 = t92+t124;
    t143 = pow(t59-t116,2.0);
    t144 = pow(t58-t117,2.0);
    t148 = -1.0/(t59-t116);
    t152 = t103+t125;
    t121 = t120*t120*t120;
    t128 = t126+1.0;
    t129 = 1.0/sqrt(t109);
    t133 = e*t25*t123;
    t134 = e*t24*t122;
    t139 = t2*t6*t19*t122;
    t140 = t2*t6*t19*t123;
    t141 = 1.0/t138;
    t145 = pow(t132,t45);
    t149 = 1.0/t143;
    t150 = t2*t6*t60*t122;
    t151 = t2*t6*t60*t123;
    t153 = pow(t132,t75);
    t161 = t110+t138;
    t162 = t113+t138;
    t169 = Winc*mu*t41*t49*t91*t96*t144;
    t170 = Winc*mu*t41*t50*t91*t96*t144;
    t173 = Winc*mu*t17*t41*t49*t91*t97*t144;
    t174 = (b_petro*e*f*sma*t11*t57*t120*-2.0)/(t59-t116);
    t180 = Wran*mu*t27*t41*t49*t96*t102*t143;
    t181 = Wran*mu*t27*t41*t50*t96*t102*t143;
    t182 = b_petro*f*t11*t57*t62*t120*t148;
    t184 = b_petro*f*t11*t57*t92*t120*t148;
    t185 = (b_petro*f*t11*t57*t62*t120)/(t59-t116);
    t186 = Wran*mu*t17*t27*t41*t49*t97*t102*t143;
    t187 = (b_petro*f*t11*t57*t92*t120)/(t59-t116);
    t146 = t13+t134;
    t147 = t10+t133;
    t159 = t76+t151;
    t160 = t77+t150;
    t163 = pow(t161,1.0/3.0);
    t164 = pow(t162,1.0/3.0);
    t171 = (Wsma*mu*t41*t51*t61*t83*t145)/4.0;
    t175 = (t141*t152)/2.0;
    t176 = (Wsma*mu*t31*t41*t51*t61*t82*t145)/4.0;
    t178 = (Wsma*mu*t31*t41*t51*t63*t83*t145)/4.0;
    t179 = Wsma*mu*t31*t41*t52*t61*t83*t145*(3.0/4.0);
    t188 = b_petro*e*f*mu*t11*t23*t57*t62*t121*t148;
    t189 = (b_petro*f*mu*sma*t11*t57*t88*t121*(-1.0/2.0))/(t59-t116);
    t190 = (b_petro*f*mu*sma*t11*t57*t88*t121)/(t59*2.0-t116*2.0);
    t206 = (Wsma*mu*n_petro*t31*t41*t45*t51*t54*t55*t61*t72*t83*t135*t153)/4.0;
    t165 = 1.0/(t163*t163);
    t166 = 1.0/(t164*t164);
    t167 = -t163;
    t168 = -t164;
    t172 = -t171;
    t177 = t31*t171;
    t183 = -t176;
    t191 = t33+t114+t175;
    t192 = Winc*mu*t41*t49*t91*t96*t147*(t58-t117)*-2.0;
    t193 = t74+t111+t175;
    t194 = Wran*mu*t27*t41*t49*t96*t102*t146*(t59-t116)*-2.0;
    t195 = b_petro*f*t11*t57*t92*t120*t146*t149;
    t196 = Winc*mu*t41*t49*t91*t96*t159*(t58-t117)*-2.0;
    t197 = Winc*mu*t41*t49*t91*t96*t159*(t58-t117)*2.0;
    t198 = Wran*mu*t27*t41*t49*t96*t102*t160*(t59-t116)*-2.0;
    t207 = -t206;
    t199 = t32+t163+t168;
    t218 = (t165*t191)/3.0;
    t220 = (t166*t193)/3.0;
    t200 = e*t199;
    t201 = t199*t199;
    t219 = -t218;
    t202 = t201-1.0;
    t203 = t200-1.0;
    t209 = t21*t88*t201;
    t210 = t23*t88*t201;
    t211 = t23*t62*t199*t200*4.0;
    t223 = t33+t219+t220;
    t204 = 1.0/t203;
    t224 = e*t223;
    t234 = sma*t21*t88*t199*t223;
    t235 = t23*t88*t199*t223*-2.0;
    t205 = t204*t204;
    t208 = t26*t204;
    t212 = t62*t204;
    t213 = t92*t204;
    t229 = t73+t164+t167+t224;
    t214 = -t213;
    t215 = t60+t212+1.0;
    t244 = -t92*t205*(t199-t224);
    t216 = t92+t214;
    t247 = t56+t208+t244;
    t217 = t216*t216;
    t228 = t202*t215*t216*2.0;
    t253 = t202*t216*(-t208+e*sma*2.0+t92*t205*(t199-t224))*-2.0;
    t221 = t202*t217;
    t242 = t209+t228;
    t245 = t199*t217*t223*2.0;
    t222 = -t221;
    t259 = t211+t235+t245+t253;
    t225 = t210+t222;
    t226 = sqrt(t225);
    t227 = 1.0/t226;
    t230 = f*mu*sma*t121*t226;
    t231 = f*t32*t120*t226;
    t232 = f*t33*t120*t226;
    t236 = (f*mu*t32*t62*t121*t226)/2.0;
    t233 = -t230;
    t237 = t187+t231;
    t254 = (f*t32*t120*t227*t242)/2.0;
    t261 = (f*t32*t120*t227*t259)/2.0;
    t238 = 1.0/(t237*t237);
    t239 = 1.0/(t237*t237*t237);
    t256 = t185+t190+t236+t254;
    t262 = -t261;
    t240 = Wape*t53*t101*t238;
    t243 = Wape*t53*t80*t84*t129*t238*2.0;
    t246 = Wape*b_petro*f*t11*t21*t53*t57*t62*t101*t120*t149*t160*t239;
    t257 = Wape*t53*t101*t239*t256*2.0;
    t265 = t174+t188+t195+t232+t233+t262;
    t266 = Wape*t53*t101*t239*(t188+t195+t232+t233+t262+b_petro*f*t11*t26*t57*t120*t148)*2.0;
    t241 = -t240;
    t255 = t197+t198+t243+t246;
    t258 = -t257;
    t267 = t157+t158+t172+t173+t178+t186+t192+t194+t266;
    t268 = t128*(t157+t158+t172+t178+t192+t194+t266+Winc*e*mu*t41*t49*t91*t97*t144*2.0+Wran*e*mu*t27*t41*t49*t97*t102*t143*2.0);
    t248 = t155+t169+t177+t180+t241;
    t260 = t156+t170+t179+t181+t183+t207+t258;
    t249 = k_petro*sma*t46*t126*t248;
    t251 = k_petro*t31*t46*t126*t248;
    t263 = t128*t260;
    t250 = -t249;
    t252 = -t251;

    D1 = t9*t92*t120*(t249-t268)-t9*t17*t23*t120*(t251-t263)+t7*t32*t92*t120*t128*t255;
    D2 = t120*(t249-t268)*(e*t124+t7*t142)-sma*t21*t69*t120*(t251-t263)+t9*t73*t120*t128*t142*t255;
    D3 = t120*t124*t128*cos(t12)*(Wape*t53*t101*t239*((b_petro*f*t92*t120*((t3/fabs(t3))))/(t59-t116)+(b_petro*f*t3*t11*t92*t120)/(t27*(t59-t116)))*-2.0+Winc*mu*t41*t49*t96*t144*(inc*2.0-inc_t*2.0)+Wran*mu*t3*t8*t41*t49*t96*t102*t143*2.0)+t3*t29*t57*t120*t124*t128*t255+Wran*mu*t8*t29*t41*t87*t94*t120*t128*t143*sin(t67)*1.0/sqrt(-t70*t70+1.0)*2.0;

    α  = atan(-D2, -D1)
    β  = atan(-D3 / sqrt(D1*D1 + D2*D2))
    return (α, β)
end