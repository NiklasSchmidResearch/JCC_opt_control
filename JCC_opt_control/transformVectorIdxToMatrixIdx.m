function matrixIdx = transformVectorIdxToMatrixIdx(matrixDims, vectorIdx)
    matrixIdx = zeros(numel(matrixDims),1);
    vectorIdx = vectorIdx - 1;
    for dim = flip(1:numel(matrixDims))
        multiplier = prod(matrixDims(1:dim-1));
        matrixIdx(dim) = 1 + floor((vectorIdx) / multiplier);
        vectorIdx = mod(vectorIdx, multiplier);
    end
end