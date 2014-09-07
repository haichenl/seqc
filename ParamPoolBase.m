classdef (Abstract) ParamPoolBase < handle
    
    methods (Access = public)
        
        % cndo2
        function SetCndo2Params(obj, cndo2params)
            if(length(cndo2params) ~= length(obj.cndo2ValidParams))
                throw(MException('ParamPoolBase:SetCndo2Params', 'Wrong number of input parameters.'));
            end
            iter = 1;
            for i = obj.cndo2ValidParams
                switch i
                    case 1
                        obj.bondingParameter = cndo2params(iter);
                    case 2
                        obj.imuAmuS = cndo2params(iter);
                    case 3
                        obj.imuAmuP = cndo2params(iter);
                    case 4
                        obj.imuAmuD = cndo2params(iter);
                    case 5
                        obj.effectiveNuclearChargeK = cndo2params(iter);
                    case 6
                        obj.effectiveNuclearChargeL = cndo2params(iter);
                    case 7
                        obj.effectiveNuclearChargeMsp = cndo2params(iter);
                    case 8
                        obj.effectiveNuclearChargeMd = cndo2params(iter);
                    case 9
                        obj.effectiveNuclearChargeNsp = cndo2params(iter);
                    otherwise
                        throw(MException('ParamPoolBase:SetCndo2Params', 'validParams wrong.'));
                end
                iter = iter + 1;
            end
        end
        function cndo2params = GetCndo2Params(obj)
            cndo2params = zeros(length(obj.cndo2ValidParams), 1);
            iter = 1;
            for i = obj.cndo2ValidParams
                switch i
                    case 1
                        cndo2params(iter) = obj.bondingParameter;
                    case 2
                        cndo2params(iter) = obj.imuAmuS;
                    case 3
                        cndo2params(iter) = obj.imuAmuP;
                    case 4
                        cndo2params(iter) = obj.imuAmuD;
                    case 5
                        cndo2params(iter) = obj.effectiveNuclearChargeK;
                    case 6
                        cndo2params(iter) = obj.effectiveNuclearChargeL;
                    case 7
                        cndo2params(iter) = obj.effectiveNuclearChargeMsp;
                    case 8
                        cndo2params(iter) = obj.effectiveNuclearChargeMd;
                    case 9
                        cndo2params(iter) = obj.effectiveNuclearChargeNsp;
                    otherwise
                        throw(MException('ParamPoolBase:GetCndo2Params', 'validParams wrong.'));
                end
                iter = iter + 1;
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
        mndoDerivedParameterD = zeros(3, 1);
        mndoDerivedParameterRho = zeros(3, 1);
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
        am1DerivedParameterD = zeros(3, 1);
        am1DerivedParameterRho = zeros(3, 1);
        am1ParameterK = zeros(4, 1);
        am1ParameterL = zeros(4, 1);
        am1ParameterM = zeros(4, 1);
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
        pm3DCoreintegralS;
        pm3DCoreintegralP;
        pm3DBondingParameterS;
        pm3DBondingParameterP;
        pm3DAlpha;
        
    end
    
end