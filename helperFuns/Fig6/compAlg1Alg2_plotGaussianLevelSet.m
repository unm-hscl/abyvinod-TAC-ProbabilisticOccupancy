function compAlg1Alg2_plotGaussianLevelSet(obstacle_mu, obstacle_sigma, ...
        level_set_threshold)
    % Plot the Gaussian level set using the contour plotting and grids

    x = -40:0.1:40;
    y = x;
    [X,Y] = meshgrid(x, y);
    XY = [X(:) Y(:)];
    p = mvnpdf(XY, obstacle_mu', obstacle_sigma); 

    contourf(X, Y, reshape(p, [], length(x)), ...
        [level_set_threshold level_set_threshold], 'FaceColor', 'k');
end
