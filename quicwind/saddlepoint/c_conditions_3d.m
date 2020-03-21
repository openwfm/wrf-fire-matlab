function [cx,cy,cz,cg]=c_conditions_3d(n,xi,yi,zi,s,cx,cy,cz,cg)
    if xi > 1
        cx((xi-1)+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(1)) = 1;
    end
    if xi < n(1)
        cx(xi+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(2)) = 1;
    end
    if yi > 1
        cy((yi-1)+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(3)) = 1;
    end
    if yi < n(2)
        cy(yi+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(4)) = 1;
    end
    if zi > 1
        cz((zi-1)+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(5)) = 1;
    end
    if zi < n(3)
        cz(zi+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(6)) = 1;
    end
    if zi==1
        cg(xi+(yi-1)*n(1),s(5)) = 1;
    end
end