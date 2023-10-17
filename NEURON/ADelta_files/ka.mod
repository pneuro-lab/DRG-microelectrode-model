COMMENT
02/09/18
Robert D. Graham
Based on channel described in Sheets et al., 2007 Journal of Physiology (https://doi.org/10.1113/jphysiol.2006.127027)
ENDCOMMENT

TITLE A-Type Potassium Channel


NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
	RANGE gkabar
	RANGE ek
	RANGE n_inf, h_inf
	RANGE tau_n, tau_h
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	gkabar = 0.07		(mho/cm2)
	ek = -68.5 			(mV)
	celsius 			(degC)
	v 					(mV)

	q10_K = 3.3
}


STATE {
	n h
}


ASSIGNED {
	ik (mA/cm2)

	tau_n
	n_inf

	tau_h
	h_inf

	q10
}


INITIAL {
	q10 = q10_K ^ ((celsius - 21) / 10)

	solveVars(v)
	n = n_inf
	h = h_inf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkabar * n * h * (v - ek)
}


DERIVATIVE states {
	solveVars(v)
	n' = (n_inf - n) / tau_n
	h' = (h_inf - h) / tau_h
}


UNITSOFF

PROCEDURE solveVars(v(mV)) {
	n_inf = (1 / (1 + exp(-(v + 5.4) / 16.4)))^4
	h_inf = 1 / (1 + exp((v + 49.9) / 4.6))
	tau_n = (0.25 + 10.04 * exp(-((v + 24.67)^2) / (2 * 34.8^2))) * (1/q10)
	tau_h = (20 + 50 * exp(-((v + 40)^2) / (2 * 40^2))) * (1/q10)

	if (tau_h < 5) {
		tau_h = 5
	}
}

UNITSON