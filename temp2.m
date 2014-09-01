% diis start
totalNumberAOs = obj.molecule.GetTotalNumberAOs();
diisStartError = Parameters.GetInstance().GetDiisStartErrorSCF();
diisEndError = Parameters.GetInstance().GetDiisEndErrorSCF();

if( 0 < diisNumErrorVect)
    diisStoredDensityMatrix(1:diisNumErrorVect-1, :) = diisStoredDensityMatrix(2:diisNumErrorVect, :);
    diisStoredDensityMatrix(diisNumErrorVect, :) = 0*diisStoredDensityMatrix(diisNumErrorVect, :);
    diisStoredErrorVect(1:diisNumErrorVect-1, :) = diisStoredErrorVect(2:diisNumErrorVect, :);
    diisStoredErrorVect(diisNumErrorVect, :) = 0*diisStoredErrorVect(diisNumErrorVect, :);
end

diisStoredDensityMatrix(diisNumErrorVect, :) = reshape(obj.orbitalElectronPopulation, 1, totalNumberAOs*totalNumberAOs);
diisStoredErrorVect(diisNumErrorVect, :) = reshape(obj.orbitalElectronPopulation - oldOrbitalElectronPopulation, 1, totalNumberAOs*totalNumberAOs);

for mi = 1:diisNumErrorVect-1
    for mj = 1:diisNumErrorVect-1
        diisErrorProducts(mi,mj) = diisErrorProducts(mi+1,mj+1);
    end
end

diisErrorProducts(diisNumErrorVect, 1:end-1) = diisStoredErrorVect * diisStoredErrorVect(diisNumErrorVect,:)';

for mi = 1:diisNumErrorVect
    diisErrorProducts(mi,diisNumErrorVect) = diisErrorProducts(diisNumErrorVect,mi);
    diisErrorProducts(mi,diisNumErrorVect+1) = -1.0;
    diisErrorProducts(diisNumErrorVect+1,mi) = -1.0;
    diisErrorCoefficients(mi) = 0.0;
end

diisErrorProducts(diisNumErrorVect+1,diisNumErrorVect+1) = 0.0;
diisErrorCoefficients(diisNumErrorVect+1) = -1.0;
diisError = max(abs(diisStoredErrorVect(diisNumErrorVect,:)));
if(diisNumErrorVect <= iterationStep && diisEndError<diisError && diisError<diisStartError)
    tmpDiisErrorProducts = diisErrorProducts;
    %                             if(abs(det(tmpDiisErrorProducts)) < 1e-12)
    %                                 hasAppliedDIIS = false;
    %                                 continue;
    %                             end
    diisErrorCoefficients = tmpDiisErrorProducts \ diisErrorCoefficients';
    obj.orbitalElectronPopulation = reshape(diisStoredDensityMatrix' * diisErrorCoefficients(1:end-1), totalNumberAOs, totalNumberAOs);
end
% diis end