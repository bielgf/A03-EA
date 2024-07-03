function plot2DBars(x,Tn,sig,units,name,section)

    % Precomputations
    X = x(:,1);
    Y = x(:,2);
    
    % Open plot window
    figure; box on; hold on; axis equal;
    % Plot deformed structure with colorbar for stresses 
    patch(X(Tn'),Y(Tn'),[sig(1,:);sig(1,:)],'EdgeColor','interp','LineWidth',2);
    % Colorbar settings
    clims = get(gca,'clim');
    clim(max(abs(clims))*[-1,1]);
    n = 128; % Number of rows
    c1 = 2/3; % Blue
    c2 = 0; % Red
    s = 0.85; % Saturation
    c = hsv2rgb([c1*ones(1,n),c2*ones(1,n);1:-1/(n-1):0,1/n:1/n:1;s*ones(1,2*n)]');
    colormap(c); 
    cb = colorbar;
    
    % Add labels
    title(sprintf('%s for %s section',name,section)); 
    xlabel('x (m)'); 
    ylabel('y (m)');
    grid on
    cb.Label.String = sprintf('Stress (%s)',units); 

end