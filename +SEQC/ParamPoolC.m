classdef (Sealed) ParamPoolC < SEQC.ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolC()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            import SEQC.Arguments;
            
            obj.coreCharge = 4.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,3,5,6];
            obj.bondingParameter = -21.0*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 14.051*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 5.572*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 5.7;
            obj.effectiveNuclearChargeL = 3.25;
            obj.effectiveNuclearChargeMsp = 0.0;
            obj.effectiveNuclearChargeMd = 0.0;
            obj.effectiveNuclearChargeNsp = 0.0;
            
            % indo
            obj.indoValidParams = [1,2,3,5,6,10,11];
            obj.indoG1 = 0.267708;
            obj.indoF2 = 0.17372;
            obj.indoF0CoefficientS = (obj.coreCharge - 0.5);
            obj.indoF0CoefficientP = (obj.coreCharge - 0.5);
            obj.indoG1CoefficientS = -1.0*(obj.coreCharge - 1.5)/6.0;
            obj.indoG1CoefficientP = -1.0/3.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = -2.0*(obj.coreCharge - 2.5)/25.0;
            
            % zindo/s
            obj.zindo_effectiveNuclearChargeMsp = 0.0;
            obj.zindo_effectiveNuclearChargeMd = 0.0;
            obj.zindoBondingParameterS = -17.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 11.11 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 55635*Arguments.GetInstance().GetKayser2AU();
            obj.zindoF2pp = 36375*Arguments.GetInstance().GetKayser2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 2;
            obj.zindoM = 2;
            obj.zindoN = 0;
            obj.zindoIonPotS = 19.84 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 10.93 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18];
            obj.mndoCoreintegralS = -52.279745 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -39.205558 * Arguments.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 1.787537;
            obj.mndoOrbitalExponentP = 1.787537;
            obj.mndoBondingParameterS = -18.985044 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -7.934122  * Arguments.GetInstance().GetEV2AU();
            obj.mndoAlpha = 2.546380 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -120.500606 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 170.89 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  12.23 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp =  11.08 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp =  11.47 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 =  9.84 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp =   2.43 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.8074661800;
            obj.mndoDerivedParameterD(3) =   0.6851577737;
            obj.mndoDerivedParameterRho(1) = 0.5/0.4494406369;
            obj.mndoDerivedParameterRho(2) = 0.5/0.6149309919;
            obj.mndoDerivedParameterRho(3) = 0.5/0.6685771472;
            
            % am1
            obj.am1ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,25,26,27,28,29,30];
            obj.am1CoreintegralS = -52.028658 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -39.614239 * Arguments.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 1.808665;
            obj.am1OrbitalExponentP = 1.685116;
            obj.am1BondingParameterS = -15.715783 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -7.719283  * Arguments.GetInstance().GetEV2AU();
            obj.am1Alpha = 2.648274 / Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss = obj.mndoGss;
            obj.am1Gpp = obj.mndoGpp;
            obj.am1Gsp = obj.mndoGsp;
            obj.am1Gpp2 = obj.mndoGpp2;
            obj.am1Hsp = obj.mndoHsp;
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.8236735591;
            obj.am1DerivedParameterD(3) = 0.7268015207;
            obj.am1DerivedParameterRho(1) = 0.5/0.4494406369;
            obj.am1DerivedParameterRho(2) = 0.5/0.6082783276;
            obj.am1DerivedParameterRho(3) = 0.5/0.6423370115;
            obj.am1ParameterK(1) = 0.011355 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.045924 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) =-0.020061 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) =-0.001260 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 5.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 1.60 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.85 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 2.05 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 2.65 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1d
            obj.am1DCoreintegralS = -52.183798 * Arguments.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = -39.368413 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -15.682341 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = -7.804762 * Arguments.GetInstance().GetEV2AU();
            obj.am1DAlpha = 2.625506 / Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3
            obj.pm3ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,27,28];
            obj.pm3CoreintegralS = -47.270320 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -36.266918 * Arguments.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 1.565085;
            obj.pm3OrbitalExponentP = 1.842345;
            obj.pm3BondingParameterS = -11.910015 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -9.802755 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Alpha = 2.707807 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.8332396384;
            obj.pm3DerivedParameterD(3) = 0.6647749859;
            obj.pm3DerivedParameterRho(1) = 0.5/0.4116151543;
            obj.pm3DerivedParameterRho(2) = 0.5/0.5885706542;
            obj.pm3DerivedParameterRho(3) = 0.5/0.7647513703;
            obj.pm3ParameterK(1) = 0.050107 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) = 0.050733 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 6.003165 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.002979 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.642214 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 0.892488 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 11.200708 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 10.796292 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gsp = 10.265027 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp2 = 9.042566 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Hsp = 2.290980 * Arguments.GetInstance().GetEV2AU();
            
            % pm3pddg
            obj.pm3PddgCoreintegralS = -48.241241 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = -36.461256 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgOrbitalExponentS = 1.567864;
            obj.pm3PddgOrbitalExponentP = 1.846659;
            obj.pm3PddgBondingParameterS = -11.952818 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP =  -9.922411 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgAlpha = 2.725772 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 0.831413;
            obj.pm3PddgDerivedParameterD(3) = 0.663222;
            obj.pm3PddgDerivedParameterRho(1) = 1.214657;
            obj.pm3PddgDerivedParameterRho(2) = 0.848467;
            obj.pm3PddgDerivedParameterRho(3) = 0.652785;
            obj.pm3PddgParameterK(1) = 0.048906 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) = 0.047697 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 5.765340 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 5.973721 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 1.682232 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 0.894406 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) =-0.000743 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) = 0.000985 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 0.836915 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 1.585236 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = -47.275431 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = -36.268916 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -11.941466 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = -9.819760 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 2.721152 / Arguments.GetInstance().GetAngstrom2AU();

        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = SEQC.ParamPoolC();
            end
            singleObj = localObj;
        end
        
    end
    
end