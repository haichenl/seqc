classdef (Sealed) ParamPoolF < ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolF()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            
            obj.coreCharge = 7.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,3,5,6];
            obj.bondingParameter = -39.0*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 32.272*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 11.080*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 8.7;
            obj.effectiveNuclearChargeL = 5.20;
            obj.effectiveNuclearChargeMsp = 0.0;
            obj.effectiveNuclearChargeMd = 0.0;
            obj.effectiveNuclearChargeNsp = 0.0;
            
            % indo
            obj.indoValidParams = [1,2,3,5,6,10,11];
            obj.indoG1 = 0.532305;
            obj.indoF2 = 0.31580;
            obj.indoF0CoefficientS = (obj.coreCharge - 0.5);
            obj.indoF0CoefficientP = (obj.coreCharge - 0.5);
            obj.indoG1CoefficientS = -1.0*(obj.coreCharge - 1.5)/6.0;
            obj.indoG1CoefficientP = -1.0/3.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = -2.0*(obj.coreCharge - 2.5)/25.0;
            
            % zindo/s
            obj.zindo_effectiveNuclearChargeMsp = 0.0;
            obj.zindo_effectiveNuclearChargeMd = 0.0;
            obj.zindoBondingParameterS = -44.0*Arguments.GetInstance().GetEV2AU(); % from orca3.0.1
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 14.00 * Arguments.GetInstance().GetEV2AU(); %  from orca3.0.1
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 116828*Arguments.GetInstance().GetKayser2AU();
            obj.zindoF2pp =  69310*Arguments.GetInstance().GetKayser2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 2;
            obj.zindoM = 5;
            obj.zindoN = 0;
            obj.zindoIonPotS = 39.39 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 20.86 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18];
            obj.mndoCoreintegralS = -131.071548 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -105.782137 * Arguments.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 2.848487;
            obj.mndoOrbitalExponentP = 2.848487;
            obj.mndoBondingParameterS = -48.290460 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -36.508540 * Arguments.GetInstance().GetEV2AU();
            obj.mndoAlpha =  3.419661 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -476.683781 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 18.86 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  16.92 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp =  16.71 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp =  17.25 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 = 14.91 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp =   4.83 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.5067166088;
            obj.mndoDerivedParameterD(3) =   0.4299633003;
            obj.mndoDerivedParameterRho(1) = 0.5/0.6217935876;
            obj.mndoDerivedParameterRho(2) = 0.5/1.0850000098;
            obj.mndoDerivedParameterRho(3) = 0.5/1.0343451703;
            
            % am1
            obj.am1ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,25,26,27,28,29,30];
            obj.am1CoreintegralS = -136.105579 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -104.889885 * Arguments.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 3.770082;
            obj.am1OrbitalExponentP = 2.494670;
            obj.am1BondingParameterS = -69.590277 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -27.922360 * Arguments.GetInstance().GetEV2AU();
            obj.am1Alpha = 5.517800/ Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss = obj.mndoGss;
            obj.am1Gpp = obj.mndoGpp;
            obj.am1Gsp = obj.mndoGsp;
            obj.am1Gpp2 = obj.mndoGpp2;
            obj.am1Hsp = obj.mndoHsp;
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.4145203025;
            obj.am1DerivedParameterD(3) = 0.4909446425;
            obj.am1DerivedParameterRho(1) = 0.5/0.6217935876;
            obj.am1DerivedParameterRho(2) = 0.5/1.2088469198;
            obj.am1DerivedParameterRho(3) = 0.5/0.9449175360;
            obj.am1ParameterK(1) = 0.242079 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.003607 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) = 0.0 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 4.80 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 4.60 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 0.930 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.660 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1DCoreintegralS = obj.am1CoreintegralS;
            obj.am1DCoreintegralP = obj.am1CoreintegralP;
            obj.am1DBondingParameterS = obj.am1BondingParameterS;
            obj.am1DBondingParameterP = obj.am1BondingParameterP;
            obj.am1DAlpha = obj.am1DAlpha;
            
            % pm3
            obj.pm3ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,27,28];
            obj.pm3CoreintegralS = -110.435303 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -105.685047 * Arguments.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 4.708555;
            obj.pm3OrbitalExponentP = 2.491178;
            obj.pm3BondingParameterS = -48.405939 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -27.744660 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Alpha = 3.358921 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.3125302275;
            obj.pm3DerivedParameterD(3) = 0.4916328225;
            obj.pm3DerivedParameterRho(1) = 0.5/0.3857423305;
            obj.pm3DerivedParameterRho(2) = 0.5/0.6768359077;
            obj.pm3DerivedParameterRho(3) = 0.5/0.6119953427;
            obj.pm3ParameterK(1) = -0.012166 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) = -0.002852 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 6.023574 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.003717 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.856859 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 2.636158 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 10.496667 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 14.817256 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gsp = 16.073689 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp2= 14.418393 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Hsp =  0.727763 * Arguments.GetInstance().GetEV2AU();
            
            % pm3pddg (unavailable yet)
            % obj.pm3PddgCoreintegralS = -48.241241 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgCoreintegralP = -36.461256 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgOrbitalExponentS = 1.567864;
            % obj.pm3PddgOrbitalExponentP = 1.846659;
            % obj.pm3PddgBondingParameterS = -11.952818 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgBondingParameterP =  -9.922411 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgAlpha = 2.725772 / Arguments.GetInstance().GetAngstrom2AU();
            % obj.pm3PddgDerivedParameterD(1) = 0.0;
            % obj.pm3PddgDerivedParameterD(2) = 0.831413;
            % obj.pm3PddgDerivedParameterD(3) = 0.663222;
            % obj.pm3PddgDerivedParameterRho(1) = 1.214657;
            % obj.pm3PddgDerivedParameterRho(2) = 0.848467;
            % obj.pm3PddgDerivedParameterRho(3) = 0.652785;
            % obj.pm3PddgParameterK(1) = 0.048906 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgParameterK(2) = 0.047697 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgParameterK(3) = 0.0;
            % obj.pm3PddgParameterK(4) = 0.0;
            % obj.pm3PddgParameterL(1) = 5.765340 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            % obj.pm3PddgParameterL(2) = 5.973721 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            % obj.pm3PddgParameterL(3) = 0.00;
            % obj.pm3PddgParameterL(4) = 0.00;
            % obj.pm3PddgParameterM(1) = 1.682232 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.pm3PddgParameterM(2) = 0.894406 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.pm3PddgParameterM(3) = 0.00;
            % obj.pm3PddgParameterM(4) = 0.00;
            % obj.pm3PddgParameterPa(1) =-0.000743 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgParameterPa(2) = 0.000985 * Arguments.GetInstance().GetEV2AU();
            % obj.pm3PddgParameterDa(1) = 0.836915 * Arguments.GetInstance().GetAngstrom2AU();
            % obj.pm3PddgParameterDa(2) = 1.585236 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = obj.pm3CoreintegralS;
            obj.pm3DCoreintegralP = obj.pm3CoreintegralP;
            obj.pm3DBondingParameterS = obj.pm3BondingParameterS;
            obj.pm3DBondingParameterP = obj.pm3BondingParameterP;
            obj.pm3DAlpha = obj.pm3Alpha;
        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = ParamPoolF();
            end
            singleObj = localObj;
        end
        
    end
    
end