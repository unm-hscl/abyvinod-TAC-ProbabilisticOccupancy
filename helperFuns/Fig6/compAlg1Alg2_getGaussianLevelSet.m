function compAlg1Alg2_getGaussianLevelSet(obstacle_mu, ...
    obstacle_sigma, obstacle_volume, alpha_value)
    % Plot the Gaussian level set using the contour plotting and grids

    x = -40:0.1:40;
    y = x;
    [X,Y] = meshgrid(x, y);
    XY = [X(:) Y(:)];
    p = mvnpdf(XY, obstacle_mu', obstacle_sigma); 

    level_set_threshold = alpha_value/obstacle_volume;
    contourf(X, Y, reshape(p, [], length(x)), ...
        [level_set_threshold level_set_threshold],'k');
end
