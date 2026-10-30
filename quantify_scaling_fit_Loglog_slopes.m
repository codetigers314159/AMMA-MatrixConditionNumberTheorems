b2 = polyfit(log(nset), log(P2(:,1)) , 1);   % slope ~ growth exponent for cond2 median
bF = polyfit(log(nset), log(PF(:,1)) , 1);   % expect ~ b2 + ~0.5-ish in practice
