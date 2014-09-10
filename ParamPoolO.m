classdef (Sealed) ParamPoolO < ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolO()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            
            obj.coreCharge = 6.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,3,5,6];
            obj.bondingParameter = -31.0*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 25.390*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 9.111*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 7.70;
            obj.effectiveNuclearChargeL = 4.55;
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
            obj.zindoBondingParameterS = -34.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 13.00 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 95298*Arguments.GetInstance().GetKayser2AU();
            obj.zindoF2pp = 55675*Arguments.GetInstance().GetKayser2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 2;
            obj.zindoM = 4;
            obj.zindoN = 0;
            obj.zindoIonPotS = 32.90 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 17.28 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoCoreintegralS = -99.64309 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -77.797472 * Arguments.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 2.699905;
            obj.mndoOrbitalExponentP = 2.699905;
            obj.mndoBondingParameterS = -32.688082 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -32.688082 * Arguments.GetInstance().GetEV2AU();
            obj.mndoAlpha = 3.160604 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -317.868506 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 59.559 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  15.42 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp =  14.52 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp =  14.48 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 = 12.98 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp =   3.94 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.5346023927;
            obj.mndoDerivedParameterD(3) =   0.4536251725;
            obj.mndoDerivedParameterRho(1) = 0.5/0.5666700426;
            obj.mndoDerivedParameterRho(2) = 0.5/0.9592303457;
            obj.mndoDerivedParameterRho(3) = 0.5/0.9495760934;
            % obj.mndoDerivedParameterD(1) =   0.0;
            % obj.mndoDerivedParameterD(2) =   0.282894 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterD(3) =   0.240043 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(1) = 0.466882 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(2) = 0.275822 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.mndoDerivedParameterRho(3) = 0.278628 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1
            obj.am1CoreintegralS = -97.830000 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -78.262380 * Arguments.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 3.108032;
            obj.am1OrbitalExponentP = 2.524039;
            obj.am1BondingParameterS = -29.272773 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -29.272773 * Arguments.GetInstance().GetEV2AU();
            obj.am1Alpha = 4.455371 / Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss = obj.mndoGss;
            obj.am1Gpp = obj.mndoGpp;
            obj.am1Gsp = obj.mndoGsp;
            obj.am1Gpp2 = obj.mndoGpp2;
            obj.am1Hsp = obj.mndoHsp;
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.4988896404;
            obj.am1DerivedParameterD(3) = 0.4852321503;
            obj.am1DerivedParameterRho(1) = 0.5/0.5666700426;
            obj.am1DerivedParameterRho(2) = 0.5/0.9960801167;
            obj.am1DerivedParameterRho(3) = 0.5/0.9065055775;
            obj.am1ParameterK(1) = 0.280962 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.081430 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) = 0.00 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) = 0.00 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 7.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 0.847918 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.445071 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1d
            obj.am1DCoreintegralS = -97.610588 * Arguments.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = -78.589700 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -29.502481 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = -29.495380 * Arguments.GetInstance().GetEV2AU();
            obj.am1DAlpha = 4.633699 / Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3
            obj.pm3CoreintegralS = -86.993002 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -71.879580 * Arguments.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 3.796544;
            obj.pm3OrbitalExponentP = 2.389402;
            obj.pm3BondingParameterS = -45.202651 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -24.752515 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Alpha = 3.217102 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.4086173087;
            obj.pm3DerivedParameterD(3) = 0.5125738036;
            obj.pm3DerivedParameterRho(1) = 0.5/0.5790088969;
            obj.pm3DerivedParameterRho(2) = 0.5/0.5299517372;
            obj.pm3DerivedParameterRho(3) = 0.5/0.8179482975;
            obj.pm3ParameterK(1) = -1.131128 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) = 1.137891 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 6.002477 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 5.950512 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.607311 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 1.598395 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 15.755760 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 13.654016 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gsp = 10.621160 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp2 = 12.40609 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Hsp = 0.593883 * Arguments.GetInstance().GetEV2AU();
            
            % pm3pddg
            obj.pm3PddgCoreintegralS = -87.412505 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = -72.183070 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgOrbitalExponentS = 3.814565;
            obj.pm3PddgOrbitalExponentP = 2.318011;
            obj.pm3PddgBondingParameterS = -44.874553 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP = -24.601939 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgAlpha = 3.225309 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 0.403741;
            obj.pm3PddgDerivedParameterD(3) = 0.528360;
            obj.pm3PddgDerivedParameterRho(1) = 0.863494;
            obj.pm3PddgDerivedParameterRho(2) = 0.936266;
            obj.pm3PddgDerivedParameterRho(3) = 0.624291;
            obj.pm3PddgParameterK(1) =-1.138455 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) = 1.146007 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 6.000043 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 5.963494 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 1.622362 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 1.614788 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) =-0.001000 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) =-0.001522 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 1.360685 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 1.366407 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = -86.960302 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = -71.926845 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -45.234302 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = -24.788037 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 3.387806 / Arguments.GetInstance().GetAngstrom2AU();
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = ParamPoolO();
            end
            singleObj = localObj;
        end
        
    end
    
end