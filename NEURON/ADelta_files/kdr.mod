COMMENT
02/09/18
Robert D. Graham
Based on channel described in Sheets et al., 2007 Journal of Physiology (2)
ENDCOMMENT

TITLE Delayed Rectifier Potassium Channel


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
	RANGE gkdrbar
	RANGE ek
	RANGE n_inf
	RANGE tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	gkdrbar = 0.018 	(mho/cm)
	ek = -68.5 			(mV)
	celsius 			(degC)
	v 					(mV)

	q10_K = 3.3
}


STATE {
	n
}


ASSIGNED {
	ik (mA/cm2)

	alpha_n
	beta_n

	tau_n
	n_inf

	q10
}


INITIAL {
	q10 = q10_K ^ ((celsius - 21) / 10)

	solveVars(v)
	n = n_inf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkdrbar * n*n*n*n  * (v - ek)
}


DERIVATIVE states {
	solveVars(v)
	n' = (n_inf - n) / tau_n
}


UNITSOFF

PROCEDURE solveVars(v(mV)) {
	n_inf = 1 / (1 + exp(-(v + 35) / 15.4))
	if (v < -31) {
		tau_n = (1000 * (0.000688 + 1/(exp((v + 75.2)/6.5) + exp((v - 131.5) / -34.8)))) * (1/q10)
	} else {
		tau_n = (0.16 + 0.8 * exp(-0.0267 * (v + 11))) * (1/q10)
	}
}

UNITSON