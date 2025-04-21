function prior = prior_up(bb,aa,x)
prior = interp1(bb,aa,x);
check = find(isnan(prior) == 1);
prior(check) = min(prior);
prior = prior/sum(prior);


end