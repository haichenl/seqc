classdef AtomCl < Atom
    
    methods
        
        function obj = AtomCl(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
        function Set_d_orbitals(obj)
            obj.valence(5) = 5;
            obj.valence(6) = 6;
            obj.valence(7) = 7;
            obj.valence(8) = 8;
            obj.valence(9) = 9;
            for i=5:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
                obj.lVec(i) = obj.realSphericalHarmonicsIndices{i}.l;
                obj.mVec(i) = obj.realSphericalHarmonicsIndices{i}.m;
            end
        end
        
        function Unset_d_orbitals(obj)
            if(length(obj.valence) > 4)
                obj.valence = obj.valence(1:4);
                obj.realSphericalHarmonicsIndices = obj.realSphericalHarmonicsIndices{1:4};
                obj.lVec = obj.lVec(1:4);
                obj.mVec = obj.mVec(1:4);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = uint8(EnumAtom.Cl);
            obj.atomicMass = 35.453*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 7.0;
            obj.numberValenceElectrons = 7;
%             obj.valenceShellType = EnumShell.mShell;
            obj.valenceShellType = 3;
%             obj.valence{1} = EnumOrbital.s;
%             obj.valence{2} = EnumOrbital.py;
%             obj.valence{3} = EnumOrbital.pz;
%             obj.valence{4} = EnumOrbital.px;
            obj.valence = (1:4)';
            obj.lVec = zeros(length(obj.valence), 1);
            obj.mVec = zeros(length(obj.valence), 1);
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
                obj.lVec(i) = obj.realSphericalHarmonicsIndices{i}.l;
                obj.mVec(i) = obj.realSphericalHarmonicsIndices{i}.m;
            end
            obj.vdWCoefficient = 8.00*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.820*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolCl.GetInstance();
            
        end
        
    end
    
end