
function  dQn(sma, e, inc, ape, ran, tru, m, ps::qLawParams)
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

    # Keep inclination and eccentricity away from zero
    if inc < 1.0e-5
        inc = 1.0e-5
    end
    if e < 1.0e-5
        e = 1.0e-5
    end

    # Start of generated code
    if (fmod(ape, 2.0*pi) > fmod(ape_t, 2.0*pi))
		g109 = fmod(ape - ape_t, 2.0*pi) >  1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : 1.0e-6;
	else
		g109 = fmod(ape - ape_t, 2.0*pi) < -1.0e-6 ? fmod(ape - ape_t, 2.0*pi) : -1.0e-6;
    end

	g110 = cos(g109);
	g2 = acos(g110);
	g3 = e*e;
	g4 = g3 - 1.0;
	g5 = 1.0 / (e*e*e);
	g6 = 1.0 / (e*e*e*e*e*e);
	g7 = g4*g4;
	g8 = g6*g7*(1.0 / 4.0);
	g9 = g8 + 1.0 / 2.7E1;
	g10 = sqrt(g9);
	g11 = 1.0 / e;
	g13 = g4*g5*(1.0 / 2.0);
	g14 = g10 + g13;
	g15 = pow(g14, 1.0 / 3.0);
	g12 = g11 + g15 - pow(g10 - g4*g5*(1.0 / 2.0), 1.0 / 3.0);
	g17 = g10 - g13;
	g18 = pow(g17, 1.0 / 3.0);
	g19 = g11 + g15 - g18;
	g24 = sma*g4;
	g25 = e*g19;
	g26 = g25 - 1.0;
	g27 = 1.0 / g26;
	g28 = sma*g4*g27;
	g16 = g24 - g28;
	g20 = cos(ape);
	g32 = mu*sma*g4;
	g21 = 1.0 / sqrt(-g32);
	g22 = b_petro + 1.0;
	g23 = g19*g19;
	g29 = g16*g16;
	g30 = sma*sma;
	g31 = g7*g23*g30;
	g33 = g23 - 1.0;
	g36 = g29*g33;
	g34 = g31 - g36;
	g35 = sqrt(g34);
	g37 = 1.0 / (e*e);
	g38 = 1.0 / sqrt(g9);
	g39 = 1.0 / (e*e*e*e*e);
	g40 = g4*g39;
	g41 = 1.0 / (e*e*e*e*e*e*e);
	g46 = g7*g41*(3.0 / 2.0);
	g42 = g40 - g46;
	g43 = g38*g42*(1.0 / 2.0);
	g44 = 1.0 / (e*e*e*e);
	g45 = 1.0 / pow(g14, 2.0 / 3.0);
	g47 = g4*g44*(3.0 / 2.0);
	g48 = 1.0 / pow(g17, 2.0 / 3.0);
	g49 = -g37 + g43 + g47;
	g50 = g48*g49*(1.0 / 3.0);
	g51 = g37 + g43 - g47;
	g53 = g45*g51*(1.0 / 3.0);
	g52 = g37 + g50 - g53;
	g54 = cos(inc);
	g55 = fabs(g54);
	g56 = sin(inc);
	g57 = 1.0 / g56;
	g58 = g20*g20;
	g65 = g3*g58;
	g59 = -g65 + 1.0;
	g60 = sqrt(g59);
	g61 = sin(ape);
	g62 = fabs(g61);
	g66 = e*g62;
	g63 = g60 - g66;
	g64 = 1.0 / g63;
	g67 = 1.0 / pow(-g32, 3.0 / 2.0);
	g68 = 1.0 / (f*f);
	g69 = sma - sma_t;
	g70 = 1.0 / sma;
	g71 = e - e_t;
	g72 = 1.0 / g4;
	g73 = inc - inc_t;
	g74 = fabs(g20);
	g75 = g61*g61;
	g78 = g3*g75;
	g76 = -g78 + 1.0;
	g77 = 1.0 / (g4*g4);
	g79 = sqrt(g76);
	g114 = e*g74;
	g80 = g79 - g114;
	g81 = g73*g73;
	g82 = 1.0 / (sma*sma*sma);
	g83 = g69*g69;
	g84 = e - 1.0;
	g85 = e + 1.0;
	g86 = fabs(g69);
	g87 = fabs(m_petro);
	g88 = 1.0 / g87;
	g89 = fabs(sma_t);
	g90 = 1.0 / g89;
	g91 = g86*g88*g90;
	g92 = pow(g91, n_petro);
	g93 = g92 + 1.0;
	g94 = 1.0 / r_petro;
	g95 = pow(g93, g94);

	#double  g100 = ran-ran_t;
	
	if (fmod(ran, 2.0*pi) > fmod(ran_t, 2.0*pi))
		g100 = fmod(ran - ran_t, 2.0*pi) >  1.0e-6 ? fmod(ran - ran_t, 2.0*pi) : 1.0e-6;
	else
		g100 = fmod(ran - ran_t, 2.0*pi) < -1.0e-6 ? fmod(ran - ran_t, 2.0*pi) : -1.0e-6;
    end

	g101 = cos(g100);
	g96 = acos(g101);
	g97 = 1.0 / sqrt(g59);
	g98 = e*g58*g97;
	g99 = g62 + g98;
	g102 = g96*g96;
	g103 = g56*g56;
	g104 = 1.0 / rpermin;
	g105 = sma*g84*g104;
	g106 = g105 + 1.0;
	g107 = k_petro*g106;
	g108 = exp(g107);
	g111 = g2*g2;
	g112 = g22*g22;
	g113 = g71*g71;
	g115 = g80*g80;
	g116 = g63*g63;
	g117 = 1.0 / g85;
	g118 = cos(tru);
	g119 = e*g118;
	g120 = g119 + 1.0;
	g121 = 1.0 / g120;
	g122 = Wp*g108;
	g123 = g122 + 1.0;
	g124 = f*g11*g21*g35;
	g128 = b_petro*f*sma*g4*g21*g55*g57*g64;
	g125 = g124 - g128;
	g126 = 1.0 / sqrt(g34);
	g127 = 1.0 / (sma*sma);
	g129 = 1.0 / (g125*g125);
	g130 = We*mu*g68*g70*g72*g113*(1.0 / 4.0);
	g131 = Winc*mu*g68*g70*g72*g81*g115;
	g132 = Wran*mu*g68*g70*g72*g102*g103*g116;
	g133 = Wsma*mu*g68*g82*g83*g84*g95*g117*(1.0 / 4.0);
	g178 = Wape*g111*g112*g129;
	g134 = g130 + g131 + g132 + g133 - g178;
	g135 = sma*g4*g121;
	g136 = g24 + g135;
	g137 = 1.0 / sqrt(g76);
	g138 = 1.0 / (g63*g63);
	g139 = 1.0 / (g125*g125*g125);

	#double  g140 = (g61 / fabs(g61));
	g140 = g61 > 0.0 ? 1.0 : -1.0;

	g141 = e*g20*g140;
	g220 = g3*g20*g61*g97;
	g142 = g141 - g220;
	g144 = sin(tru);
	g145 = f*g21*g35*g37;
	g146 = f*mu*sma*g35*g67;
	g147 = g19*g29*g52*2.0;
	g148 = e*sma*2.0;
	g149 = e*sma*g27*2.0;
	g150 = 1.0 / (g26*g26);
	g151 = e*g52;
	g152 = g11 + g15 - g18 - g151;
	g153 = sma*g4*g150*g152;
	g154 = g148 - g149 + g153;
	g155 = g16*g33*g154*2.0;
	g156 = g7*g19*g30*g52*2.0;
	g157 = e*g4*g23*g30*4.0;
	g158 = g147 - g155 - g156 + g157;
	g159 = f*g11*g21*g126*g158*(1.0 / 2.0);
	g160 = b_petro*e*f*sma*g21*g55*g57*g64*2.0;
	g161 = b_petro*f*sma*g4*g21*g55*g57*g99*g138;
	g162 = b_petro*e*f*mu*g4*g30*g55*g57*g64*g67;
	g163 = g145 - g146 - g159 + g160 + g161 + g162;
	g164 = e*2.0;
	g165 = e_t*2.0;
	g166 = g164 - g165;
	g167 = We*mu*g68*g70*g72*g166*(1.0 / 4.0);
	g168 = Wsma*mu*g68*g82*g83*g95*g117*(1.0 / 4.0);
	g169 = We*e*mu*g68*g70*g77*g113*(1.0 / 2.0);
	g170 = e*g75*g137;
	g171 = g74 + g170;
	g172 = Winc*mu*g68*g70*g72*g80*g81*g171*2.0;
	g173 = Winc*e*mu*g68*g70*g77*g81*g115*2.0;
	g174 = 1.0 / (g85*g85);
	g175 = Wsma*mu*g68*g82*g83*g84*g95*g174*(1.0 / 4.0);
	g176 = Wran*mu*g63*g68*g70*g72*g99*g102*g103*2.0;
	g177 = Wran*e*mu*g68*g70*g77*g102*g103*g116*2.0;
	g179 = Wp*k_petro*sma*g104*g108*g134;
	g180 = sma*g7*g23*2.0;
	g181 = g4*g27;
	g182 = -g3 + g181 + 1.0;
	g183 = g16*g33*g182*2.0;
	g184 = g180 + g183;
	g185 = f*g11*g21*g126*g184*(1.0 / 2.0);
	g186 = f*mu*g4*g11*g35*g67*(1.0 / 2.0);
	g187 = b_petro*f*g4*g21*g55*g57*g64;
	g188 = b_petro*f*mu*sma*g7*g55*g57*g64*g67*(1.0 / 2.0);
	g189 = g185 + g186 - g187 - g188;
	g190 = Wape*g111*g112*g139*g189*2.0;
	g191 = We*mu*g68*g72*g113*g127*(1.0 / 4.0);
	g192 = Winc*mu*g68*g72*g81*g115*g127;
	g193 = Wran*mu*g68*g72*g102*g103*g116*g127;
	g194 = sma*2.0;
	g195 = sma_t*2.0;
	g196 = g194 - g195;
	g197 = Wsma*mu*g68*g82*g84*g95*g117*g196*(1.0 / 4.0);
	g198 = 1.0 / (sma*sma*sma*sma);
	g199 = Wsma*mu*g68*g83*g84*g95*g117*g198*(3.0 / 4.0);

	#double  g200 = (g69 / fabs(g69));
	g200 = g69 > 0.0 ? 1.0 : -1.0;

	g201 = g94 - 1.0;
	g202 = pow(g93, g201);
	g203 = n_petro - 1.0;
	g204 = pow(g91, g203);
	g205 = Wsma*mu*n_petro*g68*g82*g83*g84*g88*g90*g94*g117*g200*g202*g204*(1.0 / 4.0);
	g206 = -g190 + g191 + g192 + g193 - g197 + g199 - g205;
	g207 = g123*g206;
	g208 = Wp*k_petro*g84*g104*g108*g134;
	g209 = g207 - g208;
	g210 = sin(g109);
	g211 = g110*g110;
	g212 = -g211 + 1.0;
	g213 = 1.0 / sqrt(g212);
	g214 = Wape*g2*g112*g129*g210*g213*2.0;


	#double  g215 = (g20 / fabs(g20));
	g215 = g20 > 0.0 ? 1.0 : -1.0;

	g216 = e*g61*g215;
	g217 = g3*g20*g61*g137;
	g218 = g216 - g217;
	g219 = Winc*mu*g68*g70*g72*g80*g81*g218*2.0;
	g221 = Wran*mu*g63*g68*g70*g72*g102*g103*g142*2.0;
	g222 = Wape*b_petro*f*sma*g4*g21*g55*g57*g111*g112*g138*g139*g142*2.0;
	g223 = g214 - g219 + g221 + g222;
	g143 = -g21*(g118*g136 + e*sma*g4*g121)*(g179 - g123*(-g167 - g168 + g169 + g172 + g173 + g175 + g176 + g177 + Wape*g111*g112*g163*1.0 / pow(f*g11*g21*sqrt(g31 - g29*(g12*g12 - 1.0)) - b_petro*f*sma*g4*g21*g55*g57*g64, 3.0)*2.0)) - g21*g30*g120*g209*2.0 + g11*g21*g123*g136*g144*g223;
	g224 = sma*g4*g21*g144*(g179 - g123*(-g167 - g168 + g169 + g172 + g173 + g175 + g176 + g177 + Wape*g111*g112*g139*g163*2.0)) + e*g21*g30*g144*g209*2.0 + sma*g4*g11*g21*g118*g123*g223;
	g225 = ape + tru;
	g226 = sin(g225);
	g227 = sma*g4*g21*g121*g123*cos(g225)*(Wape*g111*g112*g139*(b_petro*f*sma*g4*g21*g64*((g54 / fabs(g54))) + b_petro*f*sma*g4*g21*g54*g55*1.0 / (g56*g56)*g64)*2.0 + Winc*mu*g68*g70*g72*g115*(inc*2.0 - inc_t*2.0) + Wran*mu*g54*g56*g68*g70*g72*g102*g116*2.0) + sma*g4*g21*g54*g57*g121*g123*g223*g226 + Wran*mu*g21*g56*g68*g96*g116*g121*g123*g226*sin(g100)*1.0 / sqrt(-g101*g101 + 1.0)*2.0;
	
	Qdot_n = -g143*g143 - g224*g224 - g227*g227;

    return Qdot_n
end