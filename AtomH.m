classdef AtomH < Atom
    
    methods
        
        function obj = AtomH(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = uint8(EnumAtom.H);
            obj.atomicMass = 1.00794*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 1.0;
            obj.numberValenceElectrons = 1;
%             obj.valenceShellType = EnumShell.kShell;
            obj.valenceShellType = 1;
            obj.valence = 1;
            obj.lVec = zeros(length(obj.valence), 1);
            obj.mVec = zeros(length(obj.valence), 1);
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
                obj.lVec(i) = obj.realSphericalHarmonicsIndices{i}.l;
                obj.mVec(i) = obj.realSphericalHarmonicsIndices{i}.m;
            end
            obj.vdWCoefficient = 0.16*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.110*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolH.GetInstance();
            
        end
        
    end
    
end