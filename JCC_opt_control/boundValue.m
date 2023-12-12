function val = boundValue(val, range)
    if val<range(1)
        val=range(1);
    elseif val>range(2)
        val=range(2);
    end
end