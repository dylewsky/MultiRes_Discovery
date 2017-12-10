function xy_res = ab_to_xy(x_ab) %input has 4 rows corresponding to a, b, r, theta
    xy_res = [
    x_ab(1,:) + x_ab(2,:); ...
    x_ab(3,:) .* cos(x_ab(4,:));
    x_ab(3,:) .* sin(x_ab(4,:));
    -x_ab(1,:) + x_ab(2,:).^3
    ];
end