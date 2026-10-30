nset = 5:5:100;  T = 500;                         % 500 draws per size
P1 = []; P2 = []; PF = [];
for n = nset
    K1 = zeros(T,1); K2 = zeros(T,1); KF = zeros(T,1);
    for t=1:T
        A = randn(n,n);
        K1(t) = cond(A,1);
        K2(t) = cond(A);                          % 2-norm
        KF(t) = norm(A,'fro') * norm(inv(A),'fro'); % Frobenius conditionCLC
    end
    P1 = [P1; [prctile(K1,[50 90 99])']];
    P2 = [P2; [prctile(K2,[50 90 99])']];
    PF = [PF; [prctile(KF,[50 90 99])']];
end
% overlay medians & tails
loglog(nset, P2(:,1), 'k-',  nset, P2(:,3), 'k--');     % 2-norm median & 99th
loglog(nset, PF(:,1), 'm-',  nset, PF(:,3), 'm--');     % Frobenius median & 99th
legend('n^{1/2}','500 n^3','2-norm median','2-norm 99th','Frob median','Frob 99th','Location','northwest')
