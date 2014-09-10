classdef AtomS < Atom
    
    methods
        
        function obj = AtomS(ind)
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
            end
        end
        
        function Unset_d_orbitals(obj)
            if(length(obj.valence) > 4)
                obj.valence = obj.valence(1:4);
                obj.realSphericalHarmonicsIndices = obj.realSphericalHarmonicsIndices{1:4};
            end
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = uint8(EnumAtom.S);
            obj.atomicMass = 32.066*Arguments.GetInstance().GetGMolin2AU();
            obj.coreCharge = 6.0;
            obj.numberValenceElectrons = 6;
%             obj.valenceShellType = mShell;
            obj.valenceShellType = 3;
%             obj.valence(1) = EnumOrbital.s;
%             obj.valence(2) = EnumOrbital.py;
%             obj.valence(3) = EnumOrbital.pz;
%             obj.valence(4) = EnumOrbital.px;
            obj.valence(1) = 1;
            obj.valence(2) = 2;
            obj.valence(3) = 3;
            obj.valence(4) = 4;
%             if(Arguments.GetInstance().GetCurrentTheory() == EnumTheory.CNDO2)
%                 obj.valence(5) = EnumOrbital.dxy;
%                 obj.valence(6) = EnumOrbital.dyz;
%                 obj.valence(7) = EnumOrbital.dzz;
%                 obj.valence(8) = EnumOrbital.dzx;
%                 obj.valence(9) = EnumOrbital.dxxyy;
%             end
            for i = 1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
            end
            obj.vdWCoefficient = 10.3*Arguments.GetInstance().GetJ2AU()...
                *power(Arguments.GetInstance().GetNm2AU(),6.0)...
                /Arguments.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.870*Arguments.GetInstance().GetAngstrom2AU();
            
            obj.paramPool = ParamPoolS.GetInstance();
            
        end
        
    end
    
end