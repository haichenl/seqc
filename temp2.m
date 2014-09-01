function value = GetFockDiagElement(~, atomA, indexAtomA, mu, molecule, ...
    gammaAB, orbitalElectronPopulation, ...
    atomicElectronPopulation, ~, isGuess)
firstAOIndexA = atomA.GetFirstAOIndex();
value = atomA.GetCndo2CoreIntegral(atomA.GetValence(mu-firstAOIndexA+1), ...
    gammaAB(indexAtomA, indexAtomA), ...
    isGuess);
if(~isGuess)
    atomvect = molecule.GetAtomVect();
    temp = atomicElectronPopulation(indexAtomA) ...
        -0.5*orbitalElectronPopulation(mu, mu);
    value = value + temp*gammaAB(indexAtomA, indexAtomA);
    
    temp = 0.0;
    for BB = 1:length(molecule.GetAtomVect())
        if(BB ~= indexAtomA)
            atomBB = atomvect{BB};
            temp = temp + ( atomicElectronPopulation(BB) - atomBB.GetCoreCharge()  )...
                *gammaAB(indexAtomA, BB);
        end
    end
    value = value + temp;
end
end

function value = GetFockOffDiagElement(obj, atomA, atomB, ...
    indexAtomA, indexAtomB, mu, nu, ~, gammaAB, ...
    overlapAOs, orbitalElectronPopulation, ~, isGuess)
K = obj.GetBondingAdjustParameterK(atomA.GetValenceShellType(), atomB.GetValenceShellType());
bondParameter = 0.5*K*(atomA.GetCndo2BondingParameter() + atomB.GetCndo2BondingParameter());
value =  bondParameter*overlapAOs(mu, nu);
if(~isGuess)
    value = value - 0.5*orbitalElectronPopulation(mu, nu)*gammaAB(indexAtomA, indexAtomB);
end
end