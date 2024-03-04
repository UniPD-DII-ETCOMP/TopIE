function H = my_prepareDensityFilter(F_barh,rmin)
Nvol=size(F_barh,1);
H = zeros(Nvol,Nvol);
for id = 1:Nvol
    rVec = F_barh(id,:)-F_barh;
    dist = sqrt(sum(rVec.^2,2));  
    indDist = dist <= rmin; 
    H(id,indDist) = rmin - dist(indDist);
end
H=sparse(H);
end