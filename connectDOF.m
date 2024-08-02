function Td = connectDOF(nel,ni,Tn)
Td = zeros(nel,ni*2);

    if ni == 1
        Td = Tn;
    elseif ni == 2
        for ii = 1:nel
            for jj = 1:ni
                Td(ii,jj*2-1) = nod2dof(ni,Tn(ii,jj),1);
                Td(ii,jj*2) = nod2dof(ni,Tn(ii,jj),2);
            end
        end
    elseif ni == 3
        for ii = 1:nel
            for jj = 1:(ni-1)
                Td(ii,jj*3-2) = nod2dof(ni,Tn(ii,jj),1);
                Td(ii,jj*3-1) = nod2dof(ni,Tn(ii,jj),2);
                Td(ii,jj*3) = nod2dof(ni,Tn(ii,jj),3);
            end
        end
    end


end