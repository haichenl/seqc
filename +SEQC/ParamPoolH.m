classdef (Sealed) ParamPoolH < SEQC.ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolH()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            import SEQC.Arguments;
            
            obj.coreCharge = 1.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,5];
            obj.bondingParameter = -9.0*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 7.176*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 0.0;
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 1.2; % see P78 in J. A. Pople book
            obj.effectiveNuclearChargeL = 0.0;
            obj.effectiveNuclearChargeMsp = 0.0;
            obj.effectiveNuclearChargeMd = 0.0;
            obj.effectiveNuclearChargeNsp = 0.0;
            
            % indo
            obj.indoValidParams = [1,2,5];
            obj.indoG1 = 0.0;
            obj.indoF2 = 0.0;
            obj.indoF0CoefficientS = 0.5;
            obj.indoF0CoefficientP = 0.0;
            obj.indoG1CoefficientS = 0.0;
            obj.indoG1CoefficientP = 0.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = 0.0;
            
            % zindo/s
            obj.zindo_effectiveNuclearChargeMsp = 0.0;
            obj.zindo_effectiveNuclearChargeMd = 0.0;
            obj.zindoBondingParameterS = -12.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 12.85 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 0.0;
            obj.zindoF2pp = 0.0;
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 1;
            obj.zindoM = 0;
            obj.zindoN = 0;
            obj.zindoIonPotS = 13.06 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoValidParams = [1,3,5,7,8,16];
            obj.mndoCoreintegralS = -11.906276 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = 0.0;
            obj.mndoOrbitalExponentS = 1.331967;
            obj.mndoOrbitalExponentP = 0.0;
            obj.mndoBondingParameterS = -6.989064 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = 0.0;
            obj.mndoAlpha = 2.544134 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -11.906276 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 52.102 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss = 12.848 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) = 0.0;
            obj.mndoDerivedParameterD(2) = 0.0;
            obj.mndoDerivedParameterD(3) = 0.0;
            obj.mndoDerivedParameterRho(1) = 0.5/0.4721515374;
            obj.mndoDerivedParameterRho(2) = 0.0;
            obj.mndoDerivedParameterRho(3) = 0.0;
            
            % am1
            obj.am1ValidParams = [1,3,5,7,8,16,19,20,21,22,23,24,25,27,28,29];
            obj.am1CoreintegralS = -11.396427 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = 0.0;
            obj.am1OrbitalExponentS = 1.188078;
            obj.am1OrbitalExponentP = 0.0;
            obj.am1BondingParameterS = -6.173787 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = 0.0;
            obj.am1Alpha = 2.882324 / Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss = obj.mndoGss;
            obj.am1Gpp = obj.mndoGpp;
            obj.am1Gsp = obj.mndoGsp;
            obj.am1Gpp2 = obj.mndoGpp2;
            obj.am1Hsp = obj.mndoHsp;
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.0;
            obj.am1DerivedParameterD(3) = 0.0;
            obj.am1DerivedParameterRho(1) = 0.5/0.4721515374;
            obj.am1DerivedParameterRho(2) = 0.0;
            obj.am1DerivedParameterRho(3) = 0.0;
            obj.am1ParameterK(1) = 0.122796 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.005090 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) =-0.018336 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) = 0.000000 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 2.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 1.20 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.80 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 2.10 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1d
            obj.am1DCoreintegralS = -11.223791 * Arguments.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -6.376265 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.am1DAlpha = 3.577756 / Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3
            obj.pm3ValidParams = [1,3,5,7,8,16,19,20,21,22,23,24,27,28];
            obj.pm3CoreintegralS = -13.073321 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = 0.0;
            obj.pm3OrbitalExponentS = 0.967807;
            obj.pm3OrbitalExponentP = 0.0;
            obj.pm3BondingParameterS = -5.626512 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = 0.0;
            obj.pm3Alpha = 3.356386 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.0;
            obj.pm3DerivedParameterD(3) = 0.0;
            obj.pm3DerivedParameterRho(1) = 0.5/0.5436727936;
            obj.pm3DerivedParameterRho(2) = 0.0;
            obj.pm3DerivedParameterRho(3) = 0.0;
            obj.pm3ParameterK(1) = 1.128750 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) =-1.060329 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 5.096282 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.003788 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.537465 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 1.570189 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 14.794208 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 0.0;
            obj.pm3Gsp = 0.0;
            obj.pm3Gpp2 = 0.0;
            obj.pm3Hsp = 0.0;
            
            % pm3pddg
            obj.pm3PddgCoreintegralS = -12.893272 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = 0.0;
            obj.pm3PddgOrbitalExponentS = 0.972786;
            obj.pm3PddgOrbitalExponentP = 0.0;
            obj.pm3PddgBondingParameterS = -6.152654 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP = 0.0;
            obj.pm3PddgAlpha = 3.381686 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 0.0;
            obj.pm3PddgDerivedParameterD(3) = 0.0;
            obj.pm3PddgDerivedParameterRho(1) = 0.919616;
            obj.pm3PddgDerivedParameterRho(2) = 0.0;
            obj.pm3PddgDerivedParameterRho(3) = 0.0;
            obj.pm3PddgParameterK(1) = 1.122244 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) =-1.069737 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 4.707790 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 5.857995 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 1.547099 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 1.567893 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) = 0.057193 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) =-0.034823 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 0.663395 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 1.081901 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = -13.054076 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -5.628901 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 3.417532 / Arguments.GetInstance().GetAngstrom2AU();
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = SEQC.ParamPoolH();
            end
            singleObj = localObj;
        end
        
    end
    
end