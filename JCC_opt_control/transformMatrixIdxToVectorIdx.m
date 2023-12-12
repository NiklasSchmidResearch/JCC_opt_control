function vectorIdx = transformMatrixIdxToVectorIdx(matrixDims, matrixIdx)
    vectorIdx = 1;
    for dim = flip(1:numel(matrixDims))
        multiplier = prod(matrixDims(1:dim-1));
        vectorIdx = vectorIdx + (matrixIdx(dim)-1)*multiplier;
    end
end