%B³ad drugi i nieskonczony
function [deltaSr, deltaMx]= wskBlad(A)
    deltaSr = max(sqrt(eigs(A' * A)));
    deltaMx = max(sum(abs(A), 2));
end