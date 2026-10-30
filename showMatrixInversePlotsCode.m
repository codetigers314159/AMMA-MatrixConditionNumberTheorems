function showMatrixInversePlots(n, targetCond)
    % Generate until condition number is near targetCond
    tol = 0.5; % allow factor of 2
    found = false;
    while ~found
        A = rand(n);
        cn = cond(A);
        if cn > targetCond/2 && cn < 2*targetCond
            found = true;
        end
    end
    Ainv = inv(A);

    % Raw plots
    figure;
    subplot(2,1,1); plot(A); title(sprintf('%dx%d random A (cond≈%.2e)',n,n,cn));
    subplot(2,1,2); plot(Ainv); title('Inverse A^{-1}');

    % Row norms
    rowNormsA    = sqrt(sum(abs(A).^2,2));
    rowNormsAinv = sqrt(sum(abs(Ainv).^2,2));
    figure;
    subplot(2,1,1); plot(rowNormsA,'-o'); title('Row norms of A');
    subplot(2,1,2); plot(rowNormsAinv,'-or'); title('Row norms of A^{-1}');
end

showMatrixInversePlots(100, 1e6);