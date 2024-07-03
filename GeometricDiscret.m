function [x_prim,Tm,Tn] = GeometricDiscret(data)
    N1 = data.N1;
    N2 = data.N2;
    N3 = data.N3;

    N = [N1 N2 N3]';
    n = [N1 N1+N2 N1+N2+N3]';
    x = [
            0                   0
            0                   data.h2/2
            data.d              data.h1/2
            data.d              0
        ];
    x_mid = zeros((N1+N2+N3)+1,2);
    x_prim = zeros(2*(N1+N2+N3)+1,2);
    Tm_mid = zeros((N1+N2+N3),1);
    Tm = zeros(2*(N1+N2+N3),1);
    last = 1;
    for ii = 1:3
        xx = x(ii+1,1) - x(ii,1);
        yy = x(ii+1,2) - x(ii,2);
        rad = atan(yy/xx);
        le = sqrt(xx^2 + yy^2);
        d = le/N(ii);
        for jj = last:n(ii)
            x_mid(last,:) = x(ii,:);
            x_mid(jj+1,1) = x_mid(jj,1) + d*cos(rad);
            x_mid(jj+1,2) = x_mid(jj,2) + d*sin(rad);
            Tm_mid(jj) = ii;
        end
        last = jj+1;
    end

    x_prim(1:last,1) = flipud(x_mid(1:last,1));
    x_prim(1:last,2) = -flipud(x_mid(1:last,2));
    x_prim(last+1:end,1) = x_mid(2:end,1);
    x_prim(last+1:end,2) = x_mid(2:end,2);

    Tm(1:last-1,1) = flipud(Tm_mid(1:last-1,1));
    Tm(last:end,1) = Tm_mid;
    Tn(:,1) = 1:size(x_prim,1)-1;
    Tn(:,2) = 2:size(x_prim,1);
    Tn(end,end) = 1;
end