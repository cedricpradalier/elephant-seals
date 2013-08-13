function [f,w,jc,j] = meas_mag(common,state,obs)
    acos_prime=@(x)-1/sqrt(1-x^2);
    dot = max(-1,min(1,obs.B' * rotz(pi/2)*state.R' * obs.Bref));
    f = acos(dot);
    if nargout > 1
        w = 1./(pi/6)^2;
    end

    if nargout > 2
        jc = [];
        j = zeros(1,7);
        if (abs(f)>1e-4)
            jac = - acos_prime(dot) * cross_mat(obs.B) * state.R'*obs.Bref;
        else
            jac = zeros(3,1);
        end
        j(5:7) = jac;
    end


