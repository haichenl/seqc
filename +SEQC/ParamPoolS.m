classdef (Sealed) ParamPoolS < SEQC.ParamPoolBase
    
    methods (Access = private)
        
        function obj = ParamPoolS()
            obj.SetDefaultParameters();
        end
        
        function SetDefaultParameters(obj)
            import SEQC.Arguments;
            
            obj.coreCharge = 6.0;
            
            % cndo/2
            obj.cndo2ValidParams = [1,2,3,4,5,6,7,8];
            obj.bondingParameter = -18.150*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuS = 17.650*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuP = 6.989*Arguments.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.713*Arguments.GetInstance().GetEV2AU();
            obj.effectiveNuclearChargeK = 15.70;
            obj.effectiveNuclearChargeL = 11.85;
            obj.effectiveNuclearChargeMsp = 5.45;
            obj.effectiveNuclearChargeMd = 5.45;
            obj.effectiveNuclearChargeNsp = 0.0;
            
            % indo
            obj.indoValidParams = []; % indo's not implemented for S
            obj.indoG1 = 0.0;
            obj.indoF2 = 0.0;
            obj.indoF0CoefficientS = 0.0;
            obj.indoF0CoefficientP = 0.0;
            obj.indoG1CoefficientS = 0.0;
            obj.indoG1CoefficientP = 0.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = 0.0;
            
            % if(Arguments.GetInstance().GetCurrentTheory() == ZINDOS)
            %     obj.effectiveNuclearChargeMsp = 1.925*3.0;
            %     obj.effectiveNuclearChargeMd = 1.731*3.0;
            %
            % else
            %     obj.effectiveNuclearChargeMsp = 5.45;
            %     obj.effectiveNuclearChargeMd = 5.45;
            % end
            
            % zindo/s
            obj.zindo_effectiveNuclearChargeMsp = 1.925*3.0;
            obj.zindo_effectiveNuclearChargeMd = 1.731*3.0;
            % ORCA parameter 2.8 set
            % see "ORCA 2.8"( http:%www.thch.uni-bonn.de/tc/orca/ ).
            obj.zindoBondingParameterS = -15.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD =   0.0*Arguments.GetInstance().GetEV2AU();
            obj.zindoF0ss = 10.09 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 3.0756 * Arguments.GetInstance().GetEV2AU();
            obj.zindoF2pp = 4.5377 * Arguments.GetInstance().GetEV2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            % end (ORCA 2.8 parameter set)
            obj.zindoL = 2;
            obj.zindoM = 4;
            obj.zindoN = 0;
            obj.zindoIonPotS = 21.11 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 12.39 * Arguments.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 4.11 * Arguments.GetInstance().GetEV2AU();
            
            % mndo
            obj.mndoValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18];
            obj.mndoCoreintegralS = -72.242281 * Arguments.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -56.973207 * Arguments.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 2.312962;
            obj.mndoOrbitalExponentP = 2.009146;
            obj.mndoBondingParameterS = -10.761670 * Arguments.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -10.108433 * Arguments.GetInstance().GetEV2AU();
            obj.mndoAlpha = 2.478026 / Arguments.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -226.01239 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 66.40 * Arguments.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  12.88 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp =   9.90 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGsp =  11.26 * Arguments.GetInstance().GetEV2AU();
            obj.mndoGpp2 =  8.83 * Arguments.GetInstance().GetEV2AU();
            obj.mndoHsp =   2.26 * Arguments.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.9189935137;
            obj.mndoDerivedParameterD(3) =   0.8328513971;
            obj.mndoDerivedParameterRho(1) = 0.5/0.4733275064;
            obj.mndoDerivedParameterRho(2) = 0.5/0.5544352823;
            obj.mndoDerivedParameterRho(3) = 0.5/0.5585137839;
            
            % am1
            obj.am1ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,25,26,27,28,29,30];
            obj.am1CoreintegralS = -56.694056 * Arguments.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -48.717049 * Arguments.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 2.366515;
            obj.am1OrbitalExponentP = 1.667263;
            obj.am1BondingParameterS = -3.920566 * Arguments.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -7.905278 * Arguments.GetInstance().GetEV2AU();
            obj.am1Alpha = 2.461648 / Arguments.GetInstance().GetAngstrom2AU();
            obj.am1Gss =  11.786329 * Arguments.GetInstance().GetEV2AU();
            obj.am1Gpp =  10.039308 * Arguments.GetInstance().GetEV2AU();
            obj.am1Gsp =   8.663127 * Arguments.GetInstance().GetEV2AU();
            obj.am1Gpp2 =  7.781688 * Arguments.GetInstance().GetEV2AU();
            obj.am1Hsp =   2.532137 * Arguments.GetInstance().GetEV2AU();
            obj.am1DerivedParameterD(1) = 0.0;
            obj.am1DerivedParameterD(2) = 0.9004264562;
            obj.am1DerivedParameterD(3) = 1.0036329320;
            obj.am1DerivedParameterRho(1) = 0.5/0.4331361580;
            obj.am1DerivedParameterRho(2) = 0.5/0.5906953135;
            obj.am1DerivedParameterRho(3) = 0.5/0.6454793983;
            obj.am1ParameterK(1) =-0.509195 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) =-0.011863 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) = 0.012334 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) = 0.00 * Arguments.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 4.593691 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 5.865731 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 13.557336 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 0.00 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 0.770665 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.503313 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 2.009173 * Arguments.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 0.00 * Arguments.GetInstance().GetAngstrom2AU();
            
            % am1d
            obj.am1DCoreintegralS = -57.235044 * Arguments.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = -48.307513 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -3.311308 * Arguments.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = -7.256468 * Arguments.GetInstance().GetEV2AU();
            obj.am1DAlpha = 2.309315 / Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3
            obj.pm3ValidParams = [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,...
                19,20,21,22,23,24,27,28];
            obj.pm3CoreintegralS = -49.895371 * Arguments.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -44.392583 * Arguments.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 1.891185;
            obj.pm3OrbitalExponentP = 1.658972;
            obj.pm3BondingParameterS = -8.827465 * Arguments.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -8.091415 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Alpha = 2.269706 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 1.1214312500;
            obj.pm3DerivedParameterD(3) = 1.0086487614;
            obj.pm3DerivedParameterRho(1) = 0.5/0.3294428165;
            obj.pm3DerivedParameterRho(2) = 0.5/0.6678906502;
            obj.pm3DerivedParameterRho(3) = 0.5/0.6137333700;
            obj.pm3ParameterK(1) = -0.399191 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) = -0.054899 * Arguments.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 6.000669 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.001845 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 0.962123 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 1.579944 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 8.964667 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp = 9.968164 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gsp = 6.785936 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Gpp2 = 7.970247 * Arguments.GetInstance().GetEV2AU();
            obj.pm3Hsp = 4.041836 * Arguments.GetInstance().GetEV2AU();
            
            % pm3pddg
            obj.pm3PddgCoreintegralS = -43.906366 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = -43.461348 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgOrbitalExponentS = 1.012002;
            obj.pm3PddgOrbitalExponentP = 1.876999;
            obj.pm3PddgBondingParameterS = -2.953912 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP = -8.507779 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgAlpha = 2.539751 / Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 1.006989;
            obj.pm3PddgDerivedParameterD(3) = 0.891487;
            obj.pm3PddgDerivedParameterRho(1) = 1.517625;
            obj.pm3PddgDerivedParameterRho(2) = 0.711672;
            obj.pm3PddgDerivedParameterRho(3) = 0.754336;
            obj.pm3PddgParameterK(1) =-0.330692 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) = 0.024171 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 6.000000 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 6.000000 / power(Arguments.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 0.823837 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 2.017756 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) = 0.120434 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) =-0.002663 * Arguments.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 0.672870 * Arguments.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 2.032340 * Arguments.GetInstance().GetAngstrom2AU();
            
            % pm3d
            obj.pm3DCoreintegralS = -50.249536 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = -43.968965 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -8.397415 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = -7.594232 * Arguments.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 2.234331 / Arguments.GetInstance().GetAngstrom2AU();

        end
        
    end
    
    methods (Static)
        
        function singleObj = GetInstance()
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = SEQC.ParamPoolS();
            end
            singleObj = localObj;
        end
        
    end
    
end