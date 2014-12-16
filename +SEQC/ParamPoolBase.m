classdef (Abstract) ParamPoolBase < handle
    
    methods (Access = public)
        
        % cndo2
        function SetCndo2Params(obj, cndo2params)
            if(length(cndo2params) ~= length(obj.cndo2ValidParams))
                throw(MException('ParamPoolBase:SetCndo2Params', 'Wrong number of input parameters.'));
            end
            for i = 1:length(obj.cndo2ValidParams)
                obj.SetCndo2kthParam(obj.cndo2ValidParams(i), cndo2params(i));
            end
        end
        function cndo2params = GetCndo2Params(obj)
            cndo2params = zeros(length(obj.cndo2ValidParams), 1);
            for i = 1:length(obj.cndo2ValidParams)
                cndo2params(i) = obj.GetCndo2kthParam(obj.cndo2ValidParams(i));
            end
        end
        
        % indo
        function SetIndoParams(obj, indoparams)
            if(length(indoparams) ~= length(obj.indoValidParams))
                throw(MException('ParamPoolBase:SetIndoParams', 'Wrong number of input parameters.'));
            end
            for i = 1:length(obj.indoValidParams)
                obj.SetIndokthParam(obj.indoValidParams(i), indoparams(i));
            end
        end
        function indoparams = GetIndoParams(obj)
            indoparams = zeros(length(obj.indoValidParams), 1);
            for i = 1:length(obj.indoValidParams)
                indoparams(i) = obj.GetIndokthParam(obj.indoValidParams(i));
            end
        end
        
        % mndo
        function SetMndoParams(obj, mndoparams)
            if(length(mndoparams) ~= length(obj.mndoValidParams))
                throw(MException('ParamPoolBase:SetMndoParams', 'Wrong number of input parameters.'));
            end
            for i = 1:length(obj.mndoValidParams)
                obj.SetMndokthParam(obj.mndoValidParams(i), mndoparams(i));
            end
        end
        function mndoparams = GetMndoParams(obj)
            mndoparams = zeros(length(obj.mndoValidParams), 1);
            for i = 1:length(obj.mndoValidParams)
                mndoparams(i) = obj.GetMndokthParam(obj.mndoValidParams(i));
            end
        end
        
        % am1
        function SetAm1Params(obj, am1params)
            if(length(am1params) ~= length(obj.am1ValidParams))
                throw(MException('ParamPoolBase:SetAm1Params', 'Wrong number of input parameters.'));
            end
            for i = 1:length(obj.am1ValidParams)
                obj.SetAm1kthParam(obj.am1ValidParams(i), am1params(i));
            end
        end
        function am1params = GetAm1Params(obj)
            am1params = zeros(length(obj.am1ValidParams), 1);
            for i = 1:length(obj.am1ValidParams)
                am1params(i) = obj.GetAm1kthParam(obj.am1ValidParams(i));
            end
        end
        
        % pm3
        function SetPm3Params(obj, pm3params)
            if(length(pm3params) ~= length(obj.pm3ValidParams))
                throw(MException('ParamPoolBase:SetPm3Params', 'Wrong number of input parameters.'));
            end
            for i = 1:length(obj.pm3ValidParams)
                obj.SetPm3kthParam(obj.pm3ValidParams(i), pm3params(i));
            end
        end
        function pm3params = GetPm3Params(obj)
            pm3params = zeros(length(obj.pm3ValidParams), 1);
            for i = 1:length(obj.pm3ValidParams)
                pm3params(i) = obj.GetPm3kthParam(obj.pm3ValidParams(i));
            end
        end
        
        % zindos getters
        function res = GetZindoJss(obj)
            res = obj.zindoF0ss;
        end
        function res = GetZindoJsp(obj)
            res = obj.zindoF0ss - obj.zindoG1sp/6.0;
        end
        function res = GetZindoJsd(obj)
            res = obj.zindoF0sd - obj.zindoG2sd/10.0;
        end
        function res = GetZindoJpp(obj)
            res = obj.zindoF0ss - 2.0*obj.zindoF2pp/25.0;
        end
        function res = GetZindoJpd(obj)
            res = obj.zindoF0sd - obj.zindoG1pd/15.0 - 3.0*obj.zindoG3pd/70.0;
        end
        function res = GetZindoJdd(obj)
            res = obj.zindoF0dd - 2.0*(obj.zindoF2dd + obj.zindoF4dd)/63.0;
        end
        function res = GetZindoF0ssLower(obj)
            res = obj.zindoF0ss;
        end
        function res = GetZindoF0sdLower(obj)
            res = obj.zindoF0sd;
        end
        function res = GetZindoF0ddLower(obj)
            res = obj.zindoF0dd;
        end
        function res = GetZindoG1spLower(obj)
            res = obj.zindoG1sp/3.0;
        end
        function res = GetZindoF2ppLower(obj)
            res = obj.zindoF2pp/25.0;
        end
        function res = GetZindoG2sdLower(obj)
            res = obj.zindoG2sd/5.0;
        end
        function res = GetZindoG1pdLower(obj)
            res = obj.zindoG1pd/15.0;
        end
        function res = GetZindoF2pdLower(obj)
            res = obj.zindoF2pd/35.0;
        end
        function res = GetZindoG3pdLower(obj)
            res = obj.zindoG3pd/245.0;
        end
        function res = GetZindoF2ddLower(obj)
            res = obj.zindoF2dd/49.0;
        end
        function res = GetZindoF4ddLower(obj)
            res = obj.zindoF4dd/441.0;
        end
        
    end
    
    methods (Access = private)
        
        function SetCndo2kthParam(obj, k, param)
            switch k
                case 1
                    obj.bondingParameter = param;
                case 2
                    obj.imuAmuS = param;
                case 3
                    obj.imuAmuP = param;
                case 4
                    obj.imuAmuD = param;
                case 5
                    obj.effectiveNuclearChargeK = param;
                case 6
                    obj.effectiveNuclearChargeL = param;
                case 7
                    obj.effectiveNuclearChargeMsp = param;
                case 8
                    obj.effectiveNuclearChargeMd = param;
                case 9
                    obj.effectiveNuclearChargeNsp = param;
                otherwise
                    throw(MException('ParamPoolBase:SetCndo2kthParam', 'validParams wrong.'));
            end
        end
        function param = GetCndo2kthParam(obj, k)
            switch k
                case 1
                    param = obj.bondingParameter;
                case 2
                    param = obj.imuAmuS;
                case 3
                    param = obj.imuAmuP;
                case 4
                    param = obj.imuAmuD;
                case 5
                    param = obj.effectiveNuclearChargeK;
                case 6
                    param = obj.effectiveNuclearChargeL;
                case 7
                    param = obj.effectiveNuclearChargeMsp;
                case 8
                    param = obj.effectiveNuclearChargeMd;
                case 9
                    param = obj.effectiveNuclearChargeNsp;
                otherwise
                    throw(MException('ParamPoolBase:GetCndo2kthParam', 'validParams wrong.'));
            end
        end
        
        function SetIndokthParam(obj, k, param)
            switch k
                case 1
                    obj.bondingParameter = param;
                case 2
                    obj.imuAmuS = param;
                case 3
                    obj.imuAmuP = param;
                case 4
                    obj.imuAmuD = param;
                case 5
                    obj.effectiveNuclearChargeK = param;
                case 6
                    obj.effectiveNuclearChargeL = param;
                case 7
                    obj.effectiveNuclearChargeMsp = param;
                case 8
                    obj.effectiveNuclearChargeMd = param;
                case 9
                    obj.effectiveNuclearChargeNsp = param;
                case 10
                    obj.indoG1 = param;
                case 11
                    obj.indoF2 = param;
                otherwise
                    throw(MException('ParamPoolBase:SetIndokthParam', 'validParams wrong.'));
            end
        end
        function param = GetIndokthParam(obj, k)
            switch k
                case 1
                    param = obj.bondingParameter;
                case 2
                    param = obj.imuAmuS;
                case 3
                    param = obj.imuAmuP;
                case 4
                    param = obj.imuAmuD;
                case 5
                    param = obj.effectiveNuclearChargeK;
                case 6
                    param = obj.effectiveNuclearChargeL;
                case 7
                    param = obj.effectiveNuclearChargeMsp;
                case 8
                    param = obj.effectiveNuclearChargeMd;
                case 9
                    param = obj.effectiveNuclearChargeNsp;
                case 10
                    param = obj.indoG1;
                case 11
                    param = obj.indoF2;
                otherwise
                    throw(MException('ParamPoolBase:GetIndokthParam', 'validParams wrong.'));
            end
        end
        
        function SetMndokthParam(obj, k, param)
            switch k
                case 1
                    obj.mndoCoreintegralS = param;
                case 2
                    obj.mndoCoreintegralP = param;
                case 3
                    obj.mndoOrbitalExponentS = param;
                case 4
                    obj.mndoOrbitalExponentP = param;
                case 5
                    obj.mndoBondingParameterS = param;
                case 6
                    obj.mndoBondingParameterP = param;
                case 7
                    obj.mndoAlpha = param;
                case 8
                    obj.mndoGss = param;
                case 9
                    obj.mndoGpp = param;
                case 10
                    obj.mndoGsp = param;
                case 11
                    obj.mndoGpp2 = param;
                case 12
                    obj.mndoHsp = param;
                case 13
                    obj.mndoDerivedParameterD(1) = param;
                case 14
                    obj.mndoDerivedParameterD(2) = param;
                case 15
                    obj.mndoDerivedParameterD(3) = param;
                case 16
                    obj.mndoDerivedParameterRho(1) = param;
                case 17
                    obj.mndoDerivedParameterRho(2) = param;
                case 18
                    obj.mndoDerivedParameterRho(3) = param;
                otherwise
                    throw(MException('ParamPoolBase:SetMndokthParam', 'validParams wrong.'));
            end
        end
        function param = GetMndokthParam(obj, k)
            switch k
                case 1
                    param = obj.mndoCoreintegralS;
                case 2
                    param = obj.mndoCoreintegralP;
                case 3
                    param = obj.mndoOrbitalExponentS;
                case 4
                    param = obj.mndoOrbitalExponentP;
                case 5
                    param = obj.mndoBondingParameterS;
                case 6
                    param = obj.mndoBondingParameterP;
                case 7
                    param = obj.mndoAlpha;
                case 8
                    param = obj.mndoGss;
                case 9
                    param = obj.mndoGpp;
                case 10
                    param = obj.mndoGsp;
                case 11
                    param = obj.mndoGpp2;
                case 12
                    param = obj.mndoHsp;
                case 13
                    param = obj.mndoDerivedParameterD(1);
                case 14
                    param = obj.mndoDerivedParameterD(2);
                case 15
                    param = obj.mndoDerivedParameterD(3);
                case 16
                    param = obj.mndoDerivedParameterRho(1);
                case 17
                    param = obj.mndoDerivedParameterRho(2);
                case 18
                    param = obj.mndoDerivedParameterRho(3);
                otherwise
                    throw(MException('ParamPoolBase:GetMndokthParam', 'validParams wrong.'));
            end
        end
        
        function SetAm1kthParam(obj, k, param)
            switch k
                case 1
                    obj.am1CoreintegralS = param;
                case 2
                    obj.am1CoreintegralP = param;
                case 3
                    obj.am1OrbitalExponentS = param;
                case 4
                    obj.am1OrbitalExponentP = param;
                case 5
                    obj.am1BondingParameterS = param;
                case 6
                    obj.am1BondingParameterP = param;
                case 7
                    obj.am1Alpha = param;
                case 8
                    obj.am1Gss = param;
                case 9
                    obj.am1Gpp = param;
                case 10
                    obj.am1Gsp = param;
                case 11
                    obj.am1Gpp2 = param;
                case 12
                    obj.am1Hsp = param;
                case 13
                    obj.am1DerivedParameterD(1) = param;
                case 14
                    obj.am1DerivedParameterD(2) = param;
                case 15
                    obj.am1DerivedParameterD(3) = param;
                case 16
                    obj.am1DerivedParameterRho(1) = param;
                case 17
                    obj.am1DerivedParameterRho(2) = param;
                case 18
                    obj.am1DerivedParameterRho(3) = param;
                case 19
                    obj.am1ParameterK(1) = param;
                case 20
                    obj.am1ParameterK(2) = param;
                case 21
                    obj.am1ParameterK(3) = param;
                case 22
                    obj.am1ParameterK(4) = param;
                case 23
                    obj.am1ParameterL(1) = param;
                case 24
                    obj.am1ParameterL(2) = param;
                case 25
                    obj.am1ParameterL(3) = param;
                case 26
                    obj.am1ParameterL(4) = param;
                case 27
                    obj.am1ParameterM(1) = param;
                case 28
                    obj.am1ParameterM(2) = param;
                case 29
                    obj.am1ParameterM(3) = param;
                case 30
                    obj.am1ParameterM(4) = param;
                otherwise
                    throw(MException('ParamPoolBase:SetAm1kthParam', 'validParams wrong.'));
            end
        end
        function param = GetAm1kthParam(obj, k)
            switch k
                case 1
                    param = obj.am1CoreintegralS;
                case 2
                    param = obj.am1CoreintegralP;
                case 3
                    param = obj.am1OrbitalExponentS;
                case 4
                    param = obj.am1OrbitalExponentP;
                case 5
                    param = obj.am1BondingParameterS;
                case 6
                    param = obj.am1BondingParameterP;
                case 7
                    param = obj.am1Alpha;
                case 8
                    param = obj.am1Gss;
                case 9
                    param = obj.am1Gpp;
                case 10
                    param = obj.am1Gsp;
                case 11
                    param = obj.am1Gpp2;
                case 12
                    param = obj.am1Hsp;
                case 13
                    param = obj.am1DerivedParameterD(1);
                case 14
                    param = obj.am1DerivedParameterD(2);
                case 15
                    param = obj.am1DerivedParameterD(3);
                case 16
                    param = obj.am1DerivedParameterRho(1);
                case 17
                    param = obj.am1DerivedParameterRho(2);
                case 18
                    param = obj.am1DerivedParameterRho(3);
                case 19
                    param = obj.am1ParameterK(1);
                case 20
                    param = obj.am1ParameterK(2);
                case 21
                    param = obj.am1ParameterK(3);
                case 22
                    param = obj.am1ParameterK(4);
                case 23
                    param = obj.am1ParameterL(1);
                case 24
                    param = obj.am1ParameterL(2);
                case 25
                    param = obj.am1ParameterL(3);
                case 26
                    param = obj.am1ParameterL(4);
                case 27
                    param = obj.am1ParameterM(1);
                case 28
                    param = obj.am1ParameterM(2);
                case 29
                    param = obj.am1ParameterM(3);
                case 30
                    param = obj.am1ParameterM(4);
                otherwise
                    throw(MException('ParamPoolBase:GetAm1kthParam', 'validParams wrong.'));
            end
        end
        
        function SetPm3kthParam(obj, k, param)
            switch k
                case 1
                    obj.pm3CoreintegralS = param;
                case 2
                    obj.pm3CoreintegralP = param;
                case 3
                    obj.pm3OrbitalExponentS = param;
                case 4
                    obj.pm3OrbitalExponentP = param;
                case 5
                    obj.pm3BondingParameterS = param;
                case 6
                    obj.pm3BondingParameterP = param;
                case 7
                    obj.pm3Alpha = param;
                case 8
                    obj.pm3Gss = param;
                case 9
                    obj.pm3Gpp = param;
                case 10
                    obj.pm3Gsp = param;
                case 11
                    obj.pm3Gpp2 = param;
                case 12
                    obj.pm3Hsp = param;
                case 13
                    obj.pm3DerivedParameterD(1) = param;
                case 14
                    obj.pm3DerivedParameterD(2) = param;
                case 15
                    obj.pm3DerivedParameterD(3) = param;
                case 16
                    obj.pm3DerivedParameterRho(1) = param;
                case 17
                    obj.pm3DerivedParameterRho(2) = param;
                case 18
                    obj.pm3DerivedParameterRho(3) = param;
                case 19
                    obj.pm3ParameterK(1) = param;
                case 20
                    obj.pm3ParameterK(2) = param;
                case 21
                    obj.pm3ParameterK(3) = param;
                case 22
                    obj.pm3ParameterK(4) = param;
                case 23
                    obj.pm3ParameterL(1) = param;
                case 24
                    obj.pm3ParameterL(2) = param;
                case 25
                    obj.pm3ParameterL(3) = param;
                case 26
                    obj.pm3ParameterL(4) = param;
                case 27
                    obj.pm3ParameterM(1) = param;
                case 28
                    obj.pm3ParameterM(2) = param;
                case 29
                    obj.pm3ParameterM(3) = param;
                case 30
                    obj.pm3ParameterM(4) = param;
                otherwise
                    throw(MException('ParamPoolBase:SetPm3kthParam', 'validParams wrong.'));
            end
        end
        function param = GetPm3kthParam(obj, k)
            switch k
                case 1
                    param = obj.pm3CoreintegralS;
                case 2
                    param = obj.pm3CoreintegralP;
                case 3
                    param = obj.pm3OrbitalExponentS;
                case 4
                    param = obj.pm3OrbitalExponentP;
                case 5
                    param = obj.pm3BondingParameterS;
                case 6
                    param = obj.pm3BondingParameterP;
                case 7
                    param = obj.pm3Alpha;
                case 8
                    param = obj.pm3Gss;
                case 9
                    param = obj.pm3Gpp;
                case 10
                    param = obj.pm3Gsp;
                case 11
                    param = obj.pm3Gpp2;
                case 12
                    param = obj.pm3Hsp;
                case 13
                    param = obj.pm3DerivedParameterD(1);
                case 14
                    param = obj.pm3DerivedParameterD(2);
                case 15
                    param = obj.pm3DerivedParameterD(3);
                case 16
                    param = obj.pm3DerivedParameterRho(1);
                case 17
                    param = obj.pm3DerivedParameterRho(2);
                case 18
                    param = obj.pm3DerivedParameterRho(3);
                case 19
                    param = obj.pm3ParameterK(1);
                case 20
                    param = obj.pm3ParameterK(2);
                case 21
                    param = obj.pm3ParameterK(3);
                case 22
                    param = obj.pm3ParameterK(4);
                case 23
                    param = obj.pm3ParameterL(1);
                case 24
                    param = obj.pm3ParameterL(2);
                case 25
                    param = obj.pm3ParameterL(3);
                case 26
                    param = obj.pm3ParameterL(4);
                case 27
                    param = obj.pm3ParameterM(1);
                case 28
                    param = obj.pm3ParameterM(2);
                case 29
                    param = obj.pm3ParameterM(3);
                case 30
                    param = obj.pm3ParameterM(4);
                otherwise
                    throw(MException('ParamPoolBase:GetPm3kthParam', 'validParams wrong.'));
            end
        end
        
    end
    
    properties (SetAccess = protected)
        
        % zindo uses this
        coreCharge;
        
        % cndo/2
        cndo2ValidParams;
        bondingParameter;
        imuAmuS;
        imuAmuP;
        imuAmuD;
        effectiveNuclearChargeK;
        effectiveNuclearChargeL;
        effectiveNuclearChargeMsp;
        effectiveNuclearChargeMd;
        effectiveNuclearChargeNsp;
        
        % indo
        indoValidParams;
        indoF2;
        indoG1;
        indoF0CoefficientS;
        indoF0CoefficientP;
        indoG1CoefficientS;
        indoG1CoefficientP;
        indoF2CoefficientS;
        indoF2CoefficientP;
        
        % zindo/s
        zindo_effectiveNuclearChargeMsp;
        zindo_effectiveNuclearChargeMd;
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
        
        % mndo
        mndoValidParams;
        mndoCoreintegralS;
        mndoCoreintegralP;
        mndoOrbitalExponentS;
        mndoOrbitalExponentP;
        mndoBondingParameterS;
        mndoBondingParameterP;
        mndoAlpha;
        mndoDerivedParameterD = zeros(3, 1);
        mndoDerivedParameterRho = zeros(3, 1);
        mndoElecEnergyAtom;
        mndoHeatsFormAtom;
        mndoGss;
        mndoGpp;
        mndoGsp;
        mndoGpp2;
        mndoHsp;
        
        % am1
        am1ValidParams;
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
        am1DerivedParameterD = zeros(3, 1);
        am1DerivedParameterRho = zeros(3, 1);
        am1ParameterK = zeros(4, 1);
        am1ParameterL = zeros(4, 1);
        am1ParameterM = zeros(4, 1);
        
        % am1d
        am1DCoreintegralS;
        am1DCoreintegralP;
        am1DBondingParameterS;
        am1DBondingParameterP;
        am1DAlpha;
        
        % pm3
        pm3ValidParams;
        pm3CoreintegralS;
        pm3CoreintegralP;
        pm3OrbitalExponentS;
        pm3OrbitalExponentP;
        pm3BondingParameterS;
        pm3BondingParameterP;
        pm3Alpha;
        pm3DerivedParameterD = zeros(3, 1);
        pm3DerivedParameterRho = zeros(3, 1);
        pm3ParameterK = zeros(4, 1);
        pm3ParameterL = zeros(4, 1);
        pm3ParameterM = zeros(4, 1);
        pm3Gss;
        pm3Gpp;
        pm3Gsp;
        pm3Gpp2;
        pm3Hsp;
        
        % pm3pddg
        pm3PddgCoreintegralS;
        pm3PddgCoreintegralP;
        pm3PddgOrbitalExponentS;
        pm3PddgOrbitalExponentP;
        pm3PddgBondingParameterS;
        pm3PddgBondingParameterP;
        pm3PddgAlpha;
        pm3PddgDerivedParameterD = zeros(3, 1);
        pm3PddgDerivedParameterRho = zeros(3, 1);
        pm3PddgParameterK = zeros(4, 1);
        pm3PddgParameterL = zeros(4, 1);
        pm3PddgParameterM = zeros(4, 1);
        pm3PddgParameterPa = zeros(2, 1);
        pm3PddgParameterDa = zeros(2, 1);
        
        % pm3d
        pm3DCoreintegralS;
        pm3DCoreintegralP;
        pm3DBondingParameterS;
        pm3DBondingParameterP;
        pm3DAlpha;
        
    end
    
end