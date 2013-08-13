function v = log_mat(R)
    theta = acos((trace(R)-1)/2);
    if theta < 1e-4
        v = [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/2;
    else
        v = theta * [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/(2*sin(theta));
    end
