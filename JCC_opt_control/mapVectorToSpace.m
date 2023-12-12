function map = mapVectorToSpace(vector, mask)
    map = zeros(size(mask));
    for x = 1:size(mask,1)
        for y= 1 : size(mask,2)
            map(x,y) = vector(mask(x,y));
        end
    end
end