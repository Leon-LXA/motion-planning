function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max,axis)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 3.2.1 p constraint
    Aieq_p = [];
    bieq_p = [];
    Aieq_p_1 = [];
    Aieq_p_2 = [];
    for k = 1:n_seg
        Aieq_p_k = eye(n_order+1);
        
        Aieq_p_1 = blkdiag(Aieq_p_1,Aieq_p_k);
    end
    for k = n_seg+1 : 2*n_seg
        Aieq_p_k = -eye(n_order+1);
        Aieq_p_2 = blkdiag(Aieq_p_2,Aieq_p_k);
    end
    Aieq_p = [Aieq_p_1;Aieq_p_2];
    for n = 1:n_seg
        for m = 1:8
            bieq_p((n*8+m-8),1)=corridor_range(n,2*axis);
        end
    end
    for n = n_seg+1 : 2*n_seg
        for m = 1:8
            bieq_p((n-1)*8+m,1)=-corridor_range(n-n_seg,2*axis-1);
        end
    end
    
    
    %#####################################################
    % STEP 3.2.2 v constraint   
    %Aieq_v = [];
    Aieq_v = zeros(n_seg*7,n_seg*8);
    bieq_v = [];

    for k = 1:n_seg*7
        Aieq_v(k,k)=-1;
        Aieq_v(k,k+1)=1;
        
    end
    bieq_v = v_max * ones(7*n_seg,1);
    
    %#####################################################
    % STEP 3.2.3 a constraint   
    %Aieq_a = [];
    Aieq_a =zeros(n_seg*6,n_seg*8);
    bieq_a = [];
    for k = 1:n_seg*6
        Aieq_a(k,k) = 1;
        Aieq_a(k,k+1) = -2;
        Aieq_a(k,k+2) = 1;
        
    end
    bieq_a = a_max * ones(6*n_seg,1);
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a];
    bieq = [bieq_p; bieq_v; bieq_a];
    %Aieq = Aieq_p;
    %bieq = bieq_p;
end