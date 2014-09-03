classdef (Abstract) ParamPoolBase < handle
    
    properties (SetAccess = protected)
        
        % for now treat this as a constant
        coreCharge;
        
        % cndo/2
        imuAmuS;
        imuAmuP;
        imuAmuD;
        bondingParameter;
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
    
end