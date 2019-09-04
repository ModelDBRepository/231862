
:Comment :
: Approximation of a TREK/TWIK channel outward rectifying component assuming instantaneous activation/deactivation


NEURON	{
	SUFFIX TWIK
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik
	GLOBAL vshift
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.0013 (S/cm2) 
	vshift = 0 (mV)
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
	gk = gkbar*m^4
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
	UNITSOFF 
	mInf = 1/(1+exp(-0.032*(v-(-74+vshift))))
	UNITSON
}

