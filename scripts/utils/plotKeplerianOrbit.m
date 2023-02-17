
function plotKeplerianOrbit(kep, steps, mu)
    % Grab orbital elements
    a   = kep(1);
    e   = kep(2);
    i   = kep(3);
    ran = kep(4);
    aop = kep(5);

    % Compute states
    ths = linspace(0.0, 2*pi, steps);
    rs  = (a * (1.0 - e*e)) ./ (1.0 + e .* cos(ths));
    xs  = rs .* (cos(ths + aop) .* cos(ran) - sin(ths + aop) .* cos(i) .* sin(ran));
    ys  = rs .* (cos(ths + aop) .* sin(ran) + sin(ths + aop) .* cos(i) .* cos(ran));
    zs  = rs .* (sin(ths + aop) .* sin(i));

    plot3(xs,ys,zs, "-k")
end