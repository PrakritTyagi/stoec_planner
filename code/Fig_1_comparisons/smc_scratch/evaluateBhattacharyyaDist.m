function [ BhattDistance ] = evaluateBhattacharyyaDist( distribution_B, distribution_C )

%https://en.wikipedia.org/wiki/Bhattacharyya_distance
BC = sum(sqrt(distribution_B(:).*distribution_C(:)));
BhattDistance = -log(BC);

end

