classdef HighLevel < handle
    
    properties (SetAccess = protected)
        
        energyTotal = 0;
        energyHOMO = 0;
        energyLUMO = 0;
        energyEnvironment = 0;
        
    end
    
    methods
        
        function InitializeFromMatPsi(obj, input)
            molecule = input.molecule;
            basisSet = input.basisSet;
            fieldStrengths = input.fieldStrengths;
            
            matpsi = MatPsi(molecule.MoleculeString, basisSet);
            matpsi.fix_mol();
            matpsi.set_geom(molecule.CartesianGeometryInBohrRadius());
            
            dipole = matpsi.dipole();
            fieldMatrix = ...
                fieldStrengths(1) .* dipole(:,:,1) + ...
                fieldStrengths(2) .* dipole(:,:,2) + ...
                fieldStrengths(3) .* dipole(:,:,3);
            matpsi.RHFenv(fieldMatrix);
            obj.energyTotal = matpsi.RHF_EHF();
            energiesMO = matpsi.RHF_EMO();
            obj.energyHOMO = energiesMO(matpsi.nelec/2);
            obj.energyLUMO = energiesMO(matpsi.nelec/2+1);
            obj.energyEnvironment = ...
                reshape(2.*matpsi.RHF_D(),1,[]) * reshape(fieldMatrix,[],1);
            
        end
        
    end
    
end