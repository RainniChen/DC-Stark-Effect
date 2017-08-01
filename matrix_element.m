function result = matrix_element(n1,l1,n2,l2,h)
    
    rs = 2*n2*(n2+15);
    xs = ceil(log(rs));
    
    if ((l1 == -1) || (l2 == -1))
        X1 = radial_wav_fn2(n1, l1, h, 0);
        X2 = radial_wav_fn2(n2, l2, h, 0);
        
        m_ele = 0;
        try
            for i = 1:((2*xs/h)+1)
                m_ele = m_ele + X1(i)*X2(i)*exp(3*h*(i-1-xs/h))*h;
            end
        catch
        end
    else
        X1 = radial_wav_fn(n1, l1, h, 0);
        X2 = radial_wav_fn(n2, l2, h, 0);

        m_ele = 0;
        try
            for i = 1:(xs/h +1)
                m_ele = m_ele + X1(i)*X2(i)*exp(3*h*(i-1))*h;
            end
        catch
        end
    end
    
    result = m_ele;
    %exact = 1.5*n1*(n1^2-l1^2)^(0.5)
    
end