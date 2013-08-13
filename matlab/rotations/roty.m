
function r = roty(t)
%#eml
	ct = cos(t);
	st = sin(t);
	r =    [ct	0	st
		0	1	0
		-st	0	ct];
