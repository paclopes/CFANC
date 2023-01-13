function plot_xy_p3(x,y)
    M = size(y,2);
    if M > 1
        sorted = sort(y,2);
        plot(x,sorted(:,ceil(0.50*M)), 'color', 'black');
        hold on;
        plot(x,sorted(:,ceil(0.10*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.25*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.75*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(0.90*M)), 'color', 'black', 'LineStyle',':');
        plot(x,sorted(:,ceil(1.00*M)), 'color', [0.5 0.5 0.5], 'LineStyle',':');
        hold off;
    else
        plot(x,y, 'color', 'black')
    end

end