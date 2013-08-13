
function r = rotx(t)
%#eml
	ct = cos(t);
	st = sin(t);
	r =    [1	0	0
		0	ct	-st
		0	st	ct];
