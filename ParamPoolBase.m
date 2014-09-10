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
        
        % zindos
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
                    throw(MException('ParamPoolBase:SetIndo2kthParam', 'validParams wrong.'));
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