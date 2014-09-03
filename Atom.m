classdef Atom < handle
    
    properties (SetAccess = protected)
        
        xyz;
        atomType;
        atomicMass;
        valence;
        realSphericalHarmonicsIndices;
        valenceShellType;
        firstAOIndex;
        numberValenceElectrons;
        
        vdWCoefficient;
        vdWRadii;
        imuAmuS;
        imuAmuP;
        imuAmuD;
        bondingParameter;
        coreCharge;
        effectiveNuclearChargeK;
        effectiveNuclearChargeL;
        effectiveNuclearChargeMsp;
        effectiveNuclearChargeMd;
        effectiveNuclearChargeNsp;
        
        indoF2;
        indoG1;
        indoF0CoefficientS;
        indoF0CoefficientP;
        indoG1CoefficientS;
        indoG1CoefficientP;
        indoF2CoefficientS;
        indoF2CoefficientP;
        
        zindoBondingParameterS;
        zindoBondingParameterD;
        zindoF0ss;
        zindoF0sd;
        zindoF0dd;
        zindoG1sp;
        zindoF2pp;
        zindoG2sd;
        zindoG1pd;
        zindoF2pd;
        zindoG3pd;
        zindoF2dd;
        zindoF4dd;
        zindoL;
        zindoM;
        zindoN;
        zindoIonPotS;
        zindoIonPotP;
        zindoIonPotD;
        
        mndoCoreintegralS;         
        mndoCoreintegralP;
        mndoOrbitalExponentS;
        mndoOrbitalExponentP;
        mndoBondingParameterS;
        mndoBondingParameterP;
        mndoAlpha;
        mndoDerivedParameterD = zeros(1, 3);
        mndoDerivedParameterRho = zeros(1, 3);
        mndoElecEnergyAtom;
        mndoHeatsFormAtom;
        mndoGss;
        mndoGpp;
        mndoGsp;
        mndoGpp2;
        mndoHsp;

        am1CoreintegralS;
        am1CoreintegralP;
        am1OrbitalExponentS;
        am1OrbitalExponentP;
        am1BondingParameterS;
        am1BondingParameterP;
        am1Alpha;
        am1Gss;
        am1Gpp;
        am1Gsp;
        am1Gpp2;
        am1Hsp;
        am1DerivedParameterD = zeros(1, 3);
        am1DerivedParameterRho = zeros(1, 3);
        am1ParameterK = zeros(1, 4);
        am1ParameterL = zeros(1, 4);
        am1ParameterM = zeros(1, 4);
        am1DCoreintegralS;
        am1DCoreintegralP;
        am1DBondingParameterS;
        am1DBondingParameterP;
        am1DAlpha;
        
        pm3CoreintegralS;
        pm3CoreintegralP;
        pm3OrbitalExponentS;
        pm3OrbitalExponentP;
        pm3BondingParameterS;
        pm3BondingParameterP;
        pm3Alpha;
        pm3DerivedParameterD = zeros(1, 3);
        pm3DerivedParameterRho = zeros(1, 3);
        pm3ParameterK = zeros(1, 4);
        pm3ParameterL = zeros(1, 4);
        pm3ParameterM = zeros(1, 4);
        pm3Gss;
        pm3Gpp;
        pm3Gsp;
        pm3Gpp2;
        pm3Hsp;
        pm3PddgCoreintegralS;
        pm3PddgCoreintegralP;
        pm3PddgOrbitalExponentS;
        pm3PddgOrbitalExponentP;
        pm3PddgBondingParameterS;
        pm3PddgBondingParameterP;
        pm3PddgAlpha;
        pm3PddgDerivedParameterD = zeros(1, 3);
        pm3PddgDerivedParameterRho = zeros(1, 3);
        pm3PddgParameterK = zeros(1, 4);
        pm3PddgParameterL = zeros(1, 4);
        pm3PddgParameterM = zeros(1, 4);
        pm3PddgParameterPa = zeros(1, 2);
        pm3PddgParameterDa = zeros(1, 2);
        pm3DCoreintegralS;
        pm3DCoreintegralP;
        pm3DBondingParameterS;
        pm3DBondingParameterP;
        pm3DAlpha;
    end
    
    properties (SetAccess = private)
        index;
    end
    
    methods
        
        function obj = Atom(ind)
            obj.index = ind;
        end
        
        function SetXyz(obj, xyz_in)
            obj.xyz = xyz_in;
        end
        
        function res = GetCoreMass(obj)
            res = obj.atomicMass - obj.numberValenceElectrons;
        end
        
        function res = GetValenceSize(obj)
            res = length(obj.valence);
        end
        
        function SetCoreCharge(obj, charge)
            obj.coreCharge = charge;
        end
        
        function SetFirstAOIndex(obj, firstAOIndex_in)
            obj.firstAOIndex = firstAOIndex_in;
        end
        
        function res = GetLastAOIndex(obj)
            res = obj.firstAOIndex + length(obj.valence) - 1;
        end
        
        function res = GetEffectivePrincipalQuantumNumber(~, shellType) % ShellType shellType
%             if(shellType == 1) % EnumShell.kShell)
%                 res = 1.0;
%             elseif(shellType == 2) % EnumShell.lShell)
%                 res = 2.0;
%             elseif(shellType == 3) % EnumShell.mShell)
%                 res = 3.0;
%             elseif(shellType == 4) % EnumShell.nShell)
%                 res = 3.7;
%             else
%                 throw(MException('Atom:GetEffectivePrincipalQuantumNumber', 'Shell type wrong.'));
%             end
            res = double(shellType);
            for i = 1:length(res)
                if(res(i) == 1 || res(i) == 2 || res(i) == 3)
                    % do nothing
                elseif(res(i) == 4)
                    res(i) = 3.7;
                else
                    throw(MException('Atom:GetEffectivePrincipalQuantumNumber', 'Shell type wrong.'));
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function res = GetRealAngularPartAO(~, theta, phi, orbital) % OrbitalType orbital
            switch(orbital)
                case 1 % EnumOrbital.s
                    res = power(4.0*pi,-0.5);
                case 2 % EnumOrbital.py
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*sin(phi);
                case 3 % EnumOrbital.pz
                    res = power(3.0/(4.0*pi),0.5)*cos(theta);
                case 4 % EnumOrbital.px
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*cos(phi);
                case 5 % EnumOrbital.dxy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*sin(2.0*phi);
                case 6 % EnumOrbital.dyz
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*sin(phi);
                case 7 % EnumOrbital.dzz
                    res = power(5.0/(16.0*pi),0.5)*(3.0*power(cos(theta),2.0) - 1.0);
                case 8 % EnumOrbital.dzx
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*cos(phi);
                case 9 % EnumOrbital.dxxyy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*cos(2.0*phi);
                otherwise
                    throw(MException('Atom:GetRealAngularPartAO', 'Orbital type not possible'));
            end
        end
        
        function res = GetRadialPartAO(~, dr, orbitalExponent, shell) % ShellType shell
            principalQuantumNumber = double(shell);
            temp1 = power(2.0*orbitalExponent,principalQuantumNumber+0.5);
            temp2 = power(factorial(2*principalQuantumNumber),-0.5);
            res = temp1*temp2*power(dr,principalQuantumNumber-1)*exp(-1.0*orbitalExponent*dr);
        end
        
    end
    
end
