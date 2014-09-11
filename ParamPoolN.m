classdef (Sealed) ParamPoolN < ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolN()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            
            obj.coreCharge = 5.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,3,5,6];
            obj.bondingParameter = -25.0*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 19.316*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 7.275*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 6.7;
            obj.effectiveNuclearChargeL = 3.90;
            obj.effectiveNuclearChargeMsp = 0.0;
            obj.effectiveNuclearChargeMd = 0.0;
            obj.effectiveNuclearChargeNsp = 0.0;
            
            % indo
            obj.indoValidParams = [1,2,3,5,6,10,11];
            obj.indoG1 = 0.346029;
            obj.indoF2 = 0.219055;
            obj.indoF0CoefficientS = (obj.coreCharge - 0.5);
            obj.indoF0CoefficientP = (obj.coreCharge - 0.5);
            obj.indoG1CoefficientS = -1.0*(obj.coreCharge - 1.5)/6.0;
            obj.indoG1CoefficientP = -1.0/3.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = -2.0*(obj.coreCharge - 2.5)/25.0;
            
            % zindo/s
            obj.zindo_effectiveNuclearChargeMsp = 0.0;
            obj.zindo_effectiveNuclearChargeMd = 0.0;
            obj.zindoBondingParameterS = -26.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 12.01 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 72255*Arguments.GetInstance().GetKayser2AU();
            obj.zindoF2pp = 52100*Arguments.GetInstance().GetKayser2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 2;
            obj.zindoM = 3;
            obj.zindoN = 0;
            obj.zindoIonPotS = 25.69 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 14.05 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18];
            obj.mndoCoreintegralS = -71.932122 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -57.172319 * Arguments.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 2.255614;
            obj.mndoOrbitalExponentP = 2.255614;
            obj.mndoBondingParameterS = -20.495758 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -20.495758 * Arguments.GetInstance().GetEV2AU();
            obj.mndoAlpha = 2.861342 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -202.581201 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 113.00 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  13.59 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp =  12.98 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp =  12.66 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 = 11.59 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp =   3.14 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.6399036683;
            obj.mndoDerivedParameterD(3) =   0.5429762678;
            obj.mndoDerivedParameterRho(1) = 0.5/0.4994193177;
            obj.mndoDerivedParameterRho(2) = 0.5/0.7843433156;
            obj.mndoDerivedParameterRho(3) = 0.5/0.8126295047;
            % obj.mndoDerivedParameterD(1) =   0.0;
            % obj.mndoDerivedParameterD(2) =   0.338616 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterD(3) =   0.287325 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(1) = 0.529751 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(2) = 0.337322 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(3) = 0.324853 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1
            obj.am1ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,25,26,27,28,29,30];
            obj.am1CoreintegralS = -71.860000 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -57.167581 * Arguments.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 2.315410;
            obj.am1OrbitalExponentP = 2.157940;
            obj.am1BondingParameterS = -20.299110 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -18.238666 * Arguments.GetInstance().GetEV2AU();
            obj.am1Alpha = 2.947286 / Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss = obj.mndoGss;
            obj.am1Gpp = obj.mndoGpp;
            obj.am1Gsp = obj.mndoGsp;
            obj.am1Gpp2 = obj.mndoGpp2;
            obj.am1Hsp = obj.mndoHsp;
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.6433247425;
            obj.am1DerivedParameterD(3) = 0.5675527917;
            obj.am1DerivedParameterRho(1) = 0.5/0.4994193177;
            obj.am1DerivedParameterRho(2) = 0.5/0.7820630445;
            obj.am1DerivedParameterRho(3) = 0.5/0.7883351388;
            obj.am1ParameterK(1) = 0.025251 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.028953 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) =-0.005806 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) = 0.00 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 2.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 1.50 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 2.10 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 2.40 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1d
            obj.am1DCoreintegralS = -71.997845 * Arguments.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = -57.401718 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -20.092408 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = -18.470679 * Arguments.GetInstance().GetEV2AU();
            obj.am1DAlpha = 2.968737 / Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3
            obj.pm3CoreintegralS = -49.335672 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -47.509736 * Arguments.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 2.028094;
            obj.pm3OrbitalExponentP = 2.313728;
            obj.pm3BondingParameterS = -14.062521 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -20.043848 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Alpha = 2.830545 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.6577005762;
            obj.pm3DerivedParameterD(3) = 0.5293383109;
            obj.pm3DerivedParameterRho(1) = 0.5/0.4374893746;
            obj.pm3DerivedParameterRho(2) = 0.5/0.5030877737;
            obj.pm3DerivedParameterRho(3) = 0.5/0.7364801616;
            obj.pm3ParameterK(1) = 1.501674 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) =-1.505772 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 5.901148 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.004658 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.710740 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 1.716149 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 11.904787 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 11.754672 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gsp = 7.348565 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp2 = 10.807277 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Hsp = 1.136713 * Arguments.GetInstance().GetEV2AU();
            
            % pm3pddg
            obj.pm3PddgCoreintegralS = -49.454546 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = -47.757406 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgOrbitalExponentS = 2.035807;
            obj.pm3PddgOrbitalExponentP = 2.324327;
            obj.pm3PddgBondingParameterS = -14.117230 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP = -19.938509 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgAlpha = 2.849124 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 0.654855;
            obj.pm3PddgDerivedParameterD(3) = 0.526924;
            obj.pm3PddgDerivedParameterRho(1) = 1.142818;
            obj.pm3PddgDerivedParameterRho(2) = 0.991235;
            obj.pm3PddgDerivedParameterRho(3) = 0.676704;
            obj.pm3PddgParameterK(1) = 1.513320 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) =-1.511892 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 5.904394 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 6.030014 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 1.728376 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 1.734108 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) =-0.003160 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) = 0.012501 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 1.004172 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 1.516336 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = -49.348460 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = -47.543768 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -14.068411 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = -20.039292 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 3.060404 / Arguments.GetInstance().GetAngstrom2AU();
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = ParamPoolN();
            end
            singleObj = localObj;
        end
        
    end
    
end