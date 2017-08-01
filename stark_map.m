function stark_map(narray, h, num_points, Fdc_end)
% DC field strength in V/cm units

    Ham_0 = zeros(sum(narray)); % <-- dimension comes from m = 0
    V = zeros(sum(narray));
    
    % triplet helium quantum defect values
    del_s = 0.296670343;
    del_p = 0.06835135717;    
    
    n0 = narray(1);
    % Construct the Hamiltonian Matrix
    for index = 1:size(narray,2)
        n1 = narray(index);
        n1_s = (index-1)*n0 + 0.5*(index-1)*(index-2) + 1;
        n1_f = index*n0 + 0.5*index*(index-1);
        
        for jndex = 1:size(narray,2)
            n2 = narray(jndex);
            n2_s = (jndex-1)*n0 + 0.5*(jndex-1)*(jndex-2) + 1;
            n2_f = jndex*n0 + 0.5*jndex*(jndex-1);

            % Diagonal Componenets (unperturbed energies)
            if n1 == n2
                Ham_0(n2_s:n2_f,n2_s:n2_f) = (-0.5/n2^2)*eye(n2);
                % add quantum defects
                Ham_0(n2_s,n2_s) = (-0.5/(n2-del_s)^2);
                Ham_0(n2_s+1,n2_s+1) = (-0.5/(n2-del_p)^2);
                %Ham_0(n2_s+2,n2_s+2) = (-0.5/(n2-del_s)^2);  
            end
            
            % Off Diagonal Componenets
            for i = n1_s:n1_f
                for j = n2_s:n2_f
                    l = i-n1_s;
                    ll = j-n2_s;
                        
                    if (j-n2_s) == (i-1-n1_s)

                        if l == 0
                            mel = matrix_element(n1-del_s,l,n2,ll,h);
                        elseif ll == 0
                            mel = matrix_element(n1,l,n2-del_s,ll,h);
                        elseif l == 1
                            mel = matrix_element(n1-del_p,l,n2,ll,h);
                        elseif ll == 1
                            mel = matrix_element(n1,l,n2-del_p,ll,h);
                        elseif l == 2
                            mel = matrix_element(n1-del_d,l,n2,ll,h);
                        elseif ll == 2
                            mel = matrix_element(n1,l,n2-del_d,ll,h);
                        else
                            mel = matrix_element(n1,l,n2,ll,h);
                        end
                        V(i,j) = sqrt((l^2)/(2*l+1)/(2*l-1))*mel;
                        
                    elseif (j-n2_s) == (i+1-n1_s)
                        
                        if l == 0
                            mel = matrix_element(n1-del_s,l,n2,ll,h);
                        elseif ll == 0
                            mel = matrix_element(n1,l,n2-del_s,ll,h);
                        elseif l == 1
                            mel = matrix_element(n1-del_p,l,n2,ll,h);
                        elseif ll == 1
                            mel = matrix_element(n1,l,n2-del_p,ll,h);
                        elseif l == 2
                            mel = matrix_element(n1-del_d,l,n2,ll,h);
                        elseif ll == 2
                            mel = matrix_element(n1,l,n2-del_d,ll,h);
                        else
                            mel = matrix_element(n1,l,n2,ll,h);
                        end
                        V(i,j) = sqrt(((l+1)^2)/(2*l+3)/(2*l+1))*mel;
                    end
                    
                end
            end           
        end
    end

    Ham_0 = Ham_0*(219474.603); % unit conversion to cm^-1
    V = V*(219474.603);
    
    % Compute Eigenvalues at each Field instance
    eigens = zeros(sum(narray), num_points);
    eigens(:,1) = sort(eig(Ham_0));
    
    for index = 1:(num_points-1)
        Fdc = index*(Fdc_end/num_points);
        VE = Fdc*(1.9446907e-10)*V;
        hamiltonian = Ham_0+VE;
        eigens(:,index+1) = sort(eig(hamiltonian));
    end    
    
    % Plot
    x = linspace(0,Fdc_end, num_points);
    figure('rend','painters','pos',[0 150 750 550])
    for blah = 1:sum(narray)
        plot(x,eigens(blah,:), '.-'); hold('on')
    end
    xlabel('F_{dc} (V/cm)');
    ylabel('Energy (cm^{-1})');
end
