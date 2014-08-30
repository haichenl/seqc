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
        
        function res = GetIndex(obj)
            res = obj.index;
        end
        
        function res = GetXyz(obj)
            res = obj.xyz;
        end
        
        function SetXyz(obj, xyz_in)
            obj.xyz = xyz_in;
        end
        
        function res = GetAtomType(obj)
            res = obj.atomType;
        end
        
        function res = GetAtomicMass(obj)
            res = obj.atomicMass;
        end
        
        function res = GetCoreMass(obj)
            res = obj.atomicMass - obj.numberValenceElectrons;
        end
        
        function res = GetValenceSize(obj)
            res = length(obj.valence);
        end
        
        function res = GetValence(obj, ind)
            res = obj.valence{ind};
        end
        
        function res = GetRealSphericalHarmonicsIndex(obj, valenceIndex)
            res = obj.realSphericalHarmonicsIndices{valenceIndex};
        end
        
        function res = GetVdWCoefficient(obj)
            res = obj.vdWCoefficient;
        end
        
        function res = GetVdWRadii(obj)
            res = obj.vdWRadii;
        end
        
        function res = GetCoreCharge(obj)
            res = obj.coreCharge;
        end
        
        function SetCoreCharge(obj, charge)
            obj.coreCharge = charge;
        end
        
        function res = GetFirstAOIndex(obj)
            res = obj.firstAOIndex;
        end
        
        function SetFirstAOIndex(obj, firstAOIndex_in)
            obj.firstAOIndex = firstAOIndex_in;
        end
        
        function res = GetLastAOIndex(obj)
            res = obj.firstAOIndex + length(obj.valence) - 1;
        end
        
        function res = GetValenceShellType(obj)
            res = obj.valenceShellType;
        end
        
        function res = GetNumberValenceElectrons(obj)
            res = obj.numberValenceElectrons;
        end
        
        function res = GetAtomicBasisValue(obj, xyz_in, valenceIndex, theory)
            if(length(obj.valence) < valenceIndex)
                throw(MException('Atom:GetAtomicBasisValue', 'Given valenceIndex wrong.'));
            end
            dxyz = obj.xyz - xyz_in;
            dr = norm(dxyz);
            eulerangle = EulerAngle(dxyz);
            angularPart = obj.GetRealAngularPartAO(eulerangle.beta, eulerangle.alpha, obj.valence{valenceIndex});
            orbitalExponent = obj.GetOrbitalExponent(obj.valenceShellType, obj.valence{valenceIndex}, theory);
            radialPart = obj.GetRadialPartAO(dr, orbitalExponent, obj.valenceShellType);
            res = angularPart * radialPart;
        end
        
        function res = GetBondingParameter(obj, theory, orbital) % TheoryType theory, OrbitalType orbital
            if(nargin < 2)
                theory = TheoryType.CNDO2;
                orbital = OrbitalType.s;
            end
            if(theory == TheoryType.CNDO2 || theory == TheoryType.INDO)
                res = obj.bondingParameter;
            elseif(theory == TheoryType.ZINDOS && ( orbital == s ||...
                    orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.zindoBondingParameterS;
            elseif(theory == TheoryType.ZINDOS && ( orbital == OrbitalType.dxy ||...
                    orbital == OrbitalType.dyz ||...
                    orbital == OrbitalType.dzz ||...
                    orbital == OrbitalType.dzx ||...
                    orbital == OrbitalType.dxxyy ) )
                res = obj.zindoBondingParameterD;
            elseif(theory == TheoryType.MNDO && orbital == OrbitalType.s)
                res = obj.mndoBondingParameterS;
            elseif(theory == TheoryType.MNDO && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.mndoBondingParameterP;
            elseif(theory == TheoryType.AM1 && orbital == OrbitalType.s)
                res = obj.am1BondingParameterS;
            elseif(theory == TheoryType.AM1 && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.am1BondingParameterP;
            elseif(theory == TheoryType.AM1D && orbital == OrbitalType.s)
                res = obj.am1DBondingParameterS;
            elseif(theory == TheoryType.AM1D && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.am1DBondingParameterP;
            elseif(theory == TheoryType.PM3 && orbital == OrbitalType.s)
                res = obj.pm3BondingParameterS;
            elseif(theory == TheoryType.PM3 && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.pm3BondingParameterP;
            elseif(theory == TheoryType.PM3D && orbital == OrbitalType.s)
                res = obj.pm3DBondingParameterS;
            elseif(theory == TheoryType.PM3D && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.pm3DBondingParameterP;
            elseif(theory == TheoryType.PM3PDDG && orbital == OrbitalType.s)
                res = obj.pm3PddgBondingParameterS;
            elseif(theory == TheoryType.PM3PDDG && ( orbital == OrbitalType.px ||...
                    orbital == OrbitalType.py ||...
                    orbital == OrbitalType.pz ) )
                res = obj.pm3PddgBondingParameterP;
            else
                throw(MException('Atom:GetBondingParameter', 'Theory/Orbital type wrong.'));
            end
        end
        
        function res = GetOrbitalExponent(obj, shellType, orbitalType, theory) % ShellType shellType, OrbitalType orbitalType, TheoryType theory
            if(theory == TheoryType.CNDO2 || theory == TheoryType.INDO || theory == TheoryType.ZINDOS)
                if(shellType == ShellType.kShell && orbitalType == OrbitalType.s)
                    res = obj.effectiveNuclearChargeK/obj.GetEffectivePrincipalQuantumNumber(shellType);
                elseif(shellType == ShellType.lShell && (orbitalType == OrbitalType.s  || ...
                        orbitalType == OrbitalType.px || ...
                        orbitalType == OrbitalType.py || ...
                        orbitalType == OrbitalType.pz))
                    res = obj.effectiveNuclearChargeL/obj.GetEffectivePrincipalQuantumNumber(shellType);
                elseif(shellType == ShellType.mShell && (orbitalType == OrbitalType.s  || ...
                        orbitalType == OrbitalType.px || ...
                        orbitalType == OrbitalType.py || ...
                        orbitalType == OrbitalType.pz ))
                    res = obj.effectiveNuclearChargeMsp/obj.GetEffectivePrincipalQuantumNumber(shellType);
                elseif(shellType == ShellType.mShell && (orbitalType == OrbitalType.dxy  || ...
                        orbitalType == OrbitalType.dyz ||...
                        orbitalType == OrbitalType.dzz ||...
                        orbitalType == OrbitalType.dzx ||...
                        orbitalType == OrbitalType.dxxyy))
                    res = obj.effectiveNuclearChargeMd/obj.GetEffectivePrincipalQuantumNumber(shellType);
                elseif(shellType == ShellType.nShell && (orbitalType == OrbitalType.s  || ...
                        orbitalType == OrbitalType.px || ...
                        orbitalType == OrbitalType.py || ...
                        orbitalType == OrbitalType.pz ))
                    res = obj.effectiveNuclearChargeNsp/obj.GetEffectivePrincipalQuantumNumber(shellType);
                else
                    throw(MException('Atom:GetOrbitalExponent', 'CNDO2/INDO/ZINDOS Shell/Orbital type wrong.'));
                end
            elseif(theory == TheoryType.MNDO)
                if(orbitalType == OrbitalType.s)
                    res = obj.mndoOrbitalExponentS;
                elseif(orbitalType == OrbitalType.px ||...
                        orbitalType == OrbitalType.py ||...
                        orbitalType == OrbitalType.pz)
                    res = obj.mndoOrbitalExponentP;
                else
                    throw(MException('Atom:GetOrbitalExponent', 'MNDO Shell/Orbital type wrong.'));
                end
            elseif(theory == TheoryType.AM1 || theory == TheoryType.AM1D)
                if(orbitalType == OrbitalType.s)
                    res = obj.am1OrbitalExponentS;
                elseif(orbitalType == OrbitalType.px ||...
                        orbitalType == OrbitalType.py ||...
                        orbitalType == OrbitalType.pz)
                    res = obj.am1OrbitalExponentP;
                else
                    throw(MException('Atom:GetOrbitalExponent', 'AM1/AM1D Orbital type wrong.'));
                end
            elseif(theory == TheoryType.PM3 || theory == TheoryType.PM3D)
                if(orbitalType == OrbitalType.s)
                    res = obj.pm3OrbitalExponentS;
                elseif(orbitalType == OrbitalType.px ||...
                        orbitalType == OrbitalType.py ||...
                        orbitalType == OrbitalType.pz)
                    res = obj.pm3OrbitalExponentP;
                else
                    throw(MException('Atom:GetOrbitalExponent', 'PM3/PM3D Orbital type wrong.'));
                end
            elseif(theory == TheoryType.PM3PDDG)
                if(orbitalType == OrbitalType.s)
                    res = obj.pm3PddgOrbitalExponentS;
                elseif(orbitalType == OrbitalType.px ||...
                        orbitalType == OrbitalType.py ||...
                        orbitalType == OrbitalType.pz)
                    res = obj.pm3PddgOrbitalExponentP;
                else
                    throw(MException('Atom:GetOrbitalExponent', 'PM3DDG Orbital type wrong.'));
                end
            else
                throw(MException('Atom:GetOrbitalExponent', 'Theory type wrong.'));
            end
        end
        
        % todo: complete other theories
        function value = GetCoreIntegral(obj, orbital, gamma, isGuess, theory)
            if(theory == TheoryType.CNDO2)
                value = obj.GetCndo2CoreIntegral(orbital, gamma, isGuess);
            else
                throw(MException('Atom:GetCoreIntegral', 'Not implemented for this theory yet.'));
            end
        end
        
    end
    
    methods (Access = private)
        
        function res = GetRealAngularPartAO(~, theta, phi, orbital) % OrbitalType orbital
            switch(orbital)
                case OrbitalType.s
                    res = power(4.0*pi,-0.5);
                case OrbitalType.py
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*sin(phi);
                case OrbitalType.pz
                    res = power(3.0/(4.0*pi),0.5)*cos(theta);
                case OrbitalType.px
                    res = power(3.0/(4.0*pi),0.5)*sin(theta)*cos(phi);
                case OrbitalType.dxy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*sin(2.0*phi);
                case OrbitalType.dyz
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*sin(phi);
                case OrbitalType.dzz
                    res = power(5.0/(16.0*pi),0.5)*(3.0*power(cos(theta),2.0) - 1.0);
                case OrbitalType.dzx
                    res = power(15.0/(16.0*pi),0.5)*sin(2.0*theta)*cos(phi);
                case OrbitalType.dxxyy
                    res = power(15.0/(16.0*pi),0.5)*power(sin(theta),2.0)*cos(2.0*phi);
                otherwise
                    throw(MException('Atom:GetRealAngularPartAO', 'Orbital type not possible'));
            end
        end
        
        function res = GetEffectivePrincipalQuantumNumber(~, shellType) % ShellType shellType
            if(shellType == ShellType.kShell)
                res = 1.0;
            elseif(shellType == ShellType.lShell)
                res = 2.0;
            elseif(shellType == ShellType.mShell)
                res = 3.0;
            elseif(shellType == ShellType.nShell)
                res = 3.7;
            else
                throw(MException('Atom:GetEffectivePrincipalQuantumNumber', 'Shell type wrong.'));
            end
        end
        
        function res = GetRadialPartAO(~, dr, orbitalExponent, shell) % ShellType shell
            principalQuantumNumber = double(shell);
            temp1 = power(2.0*orbitalExponent,principalQuantumNumber+0.5);
            temp2 = power(factorial(2*principalQuantumNumber),-0.5);
            res = temp1*temp2*power(dr,principalQuantumNumber-1)*exp(-1.0*orbitalExponent*dr);
        end
        
        function res = GetCndo2CoreIntegral(obj, orbital, gamma, isGuess) % OrbitalType orbital
            if(orbital == OrbitalType.s)
                res = -1.0*obj.imuAmuS;
            elseif(orbital == OrbitalType.px || orbital == OrbitalType.py || orbital == OrbitalType.pz)
                res = -1.0*obj.imuAmuP;
            elseif(orbital == OrbitalType.dxy || ...
                    orbital == OrbitalType.dyz || ...
                    orbital == OrbitalType.dzz || ...
                    orbital == OrbitalType.dzx || ...
                    orbital == OrbitalType.dxxyy )
                res = -1.0*obj.imuAmuD;
            else
                throw(MException('Atom:GetCndo2CoreIntegral', 'CNDO2 Orbital type wrong.'));
            end
            if(~isGuess)
                res = res - (obj.coreCharge - 0.5)*gamma;
            end
        end
        
    end
    
end
