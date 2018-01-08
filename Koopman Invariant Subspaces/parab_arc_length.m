function xFlat = parab_arc_length(x)
    xs = sign(x);
    x = abs(x);
    xFlatAbs = 0.5 * x .* (1 + 4*x.^2).^(1/2) + 0.25 * log(2*x + (1 + 4*x.^2).^(1/2));
    xFlat = xFlatAbs .* xs;
end