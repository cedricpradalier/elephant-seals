function Sinv = Srot_inv(R)
    c1 = R(:,1); c2 = R(:,2); c3 = R(:,3); 
    c1_x = cross_mat(c1);
    c2_x = cross_mat(c2);
    c3_x = cross_mat(c3);
    Sinv = [c1_x c2_x c3_x];

