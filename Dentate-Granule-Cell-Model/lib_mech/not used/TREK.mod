
:Comment :
: Approximation of a TREK/TWIK channel outward rectifying component assuming instantaneous activation/deactivation


NEURON	{
	SUFFIX TREK
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik
	GLOBAL vshift, b, d
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.0013 (S/cm2) 
	vshift = 0 (mV)
	b = 0.1 (/mV)
	d = 0 (mV)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	mInf
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gk = gkbar*m
	ik = gk*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
	LOCAL x
	UNITSOFF 
	x = b*(v-d)
	if (fabs(x) > 1e-6) {
		mInf = x/(1-exp(-x))
	}else{
		mInf = 1/(1-0.5*x)
	}
	UNITSON
}

