function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        M_k = [];
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        %
        %
        %
        %
        M_k(1:4,1:8)= [1 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0;0 0 2 0 0 0 0 0; 0 0 0 6 0 0 0 0];
        %M_k(5,1:8)=tvec(ts(k),n_order,0);
        for r = 0:3
            M_k(r+5,1:8) = calc_tvec(ts(k),n_order,r);
        end
        M = blkdiag(M, M_k);
    end
end