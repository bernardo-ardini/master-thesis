function out=asmbvec(int,k,j,i,space)
    va=int.*reshape(space.geo.gauss{k}.trace{j}.measure,1,[]);
    va=va(:);
    ro=space.ind{i};
    ro=ro(:);
    mask=ro>0;
    out=accumarray(ro(mask),va(mask),[space.ndof,1]);
end