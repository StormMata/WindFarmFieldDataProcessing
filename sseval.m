function sse = sseval(x,xdata,ydata)
%     A = x(1);
%     lambda = x(2);
%     sse = sum((ydata - A*exp(-lambda*tdata)).^2);

    uref = x(1);
    a = x(2);

    sse = sum((ydata - uref.*(xdata/43).^a).^2);
end