classdef CAtom < Atom
    
    methods
        
        function obj = CAtom(ind)
            obj@Atom(ind);
            obj.SetAtomicParameters();
        end
        
    end
    
    methods (Access = protected)
        
        function SetAtomicParameters(obj)
            obj.atomType = AtomType.C;
            obj.atomicMass = 12.0107*Parameters.GetInstance().GetGMolin2AU();
            obj.coreCharge = 4.0;
            obj.numberValenceElectrons = 4;
%             obj.valenceShellType = ShellType.lShell;
            obj.valenceShellType = 2;
%             obj.valence{1} = OrbitalType.s;
%             obj.valence{2} = OrbitalType.py;
%             obj.valence{3} = OrbitalType.pz;
%             obj.valence{4} = OrbitalType.px;
            obj.valence(1) = 1;
            obj.valence(2) = 2;
            obj.valence(3) = 3;
            obj.valence(4) = 4;
            for i=1:length(obj.valence)
                obj.realSphericalHarmonicsIndices{i} = RealSphericalHarmonicsIndex(obj.valence(i));
            end
            obj.vdWCoefficient = 1.65*Parameters.GetInstance().GetJ2AU()...
                *power(Parameters.GetInstance().GetNm2AU(),6.0)...
                /Parameters.GetInstance().GetAvogadro();
            obj.vdWRadii = 1.610*Parameters.GetInstance().GetAngstrom2AU();
            obj.bondingParameter = -21.0*Parameters.GetInstance().GetEV2AU();
            obj.imuAmuS = 14.051*Parameters.GetInstance().GetEV2AU();
            obj.imuAmuP = 5.572*Parameters.GetInstance().GetEV2AU();
            obj.imuAmuD = 0.0;
            obj.effectiveNuclearChargeK = 5.7;
            obj.effectiveNuclearChargeL = 3.25;
            obj.effectiveNuclearChargeMsp = 0.0;
            obj.effectiveNuclearChargeMd = 0.0;
            obj.indoG1 = 0.267708;
            obj.indoF2 = 0.17372;
            obj.indoF0CoefficientS = (obj.coreCharge - 0.5);
            obj.indoF0CoefficientP = (obj.coreCharge - 0.5);
            obj.indoG1CoefficientS = -1.0*(obj.coreCharge - 1.5)/6.0;
            obj.indoG1CoefficientP = -1.0/3.0;
            obj.indoF2CoefficientS = 0.0;
            obj.indoF2CoefficientP = -2.0*(obj.coreCharge - 2.5)/25.0;
            obj.zindoBondingParameterS = -17.0*Parameters.GetInstance().GetEV2AU();
            obj.zindoBondingParameterD = 0.0;
            obj.zindoF0ss = 11.11 * Parameters.GetInstance().GetEV2AU();
            obj.zindoF0sd = 0.0;
            obj.zindoF0dd = 0.0;
            obj.zindoG1sp = 55635*Parameters.GetInstance().GetKayser2AU();
            obj.zindoF2pp = 36375*Parameters.GetInstance().GetKayser2AU();
            obj.zindoG2sd = 0.0;
            obj.zindoG1pd = 0.0;
            obj.zindoF2pd = 0.0;
            obj.zindoG3pd = 0.0;
            obj.zindoF2dd = 0.0;
            obj.zindoF4dd = 0.0;
            obj.zindoL = 2;
            obj.zindoM = 2;
            obj.zindoN = 0;
            obj.zindoIonPotS = 19.84 * Parameters.GetInstance().GetEV2AU();
            obj.zindoIonPotP = 10.93 * Parameters.GetInstance().GetEV2AU();
            obj.zindoIonPotD = 0.0 * Parameters.GetInstance().GetEV2AU();
            obj.mndoCoreintegralS = -52.279745 * Parameters.GetInstance().GetEV2AU();
            obj.mndoCoreintegralP = -39.205558 * Parameters.GetInstance().GetEV2AU();
            obj.mndoOrbitalExponentS = 1.787537;
            obj.mndoOrbitalExponentP = 1.787537;
            obj.mndoBondingParameterS = -18.985044 * Parameters.GetInstance().GetEV2AU();
            obj.mndoBondingParameterP = -7.934122  * Parameters.GetInstance().GetEV2AU();
            obj.mndoAlpha = 2.546380 / Parameters.GetInstance().GetAngstrom2AU();
            obj.mndoElecEnergyAtom = -120.500606 * Parameters.GetInstance().GetEV2AU();
            obj.mndoHeatsFormAtom = 170.89 * Parameters.GetInstance().GetKcalMolin2AU();
            obj.mndoGss =  12.23 * Parameters.GetInstance().GetEV2AU();
            obj.mndoGpp =  11.08 * Parameters.GetInstance().GetEV2AU();
            obj.mndoGsp =  11.47 * Parameters.GetInstance().GetEV2AU();
            obj.mndoGpp2 =  9.84 * Parameters.GetInstance().GetEV2AU();
            obj.mndoHsp =   2.43 * Parameters.GetInstance().GetEV2AU();
            obj.mndoDerivedParameterD(1) =   0.0;
            obj.mndoDerivedParameterD(2) =   0.8074661800;
            obj.mndoDerivedParameterD(3) =   0.6851577737;
            obj.mndoDerivedParameterRho(1) = 0.5/0.4494406369;
            obj.mndoDerivedParameterRho(2) = 0.5/0.6149309919;
            obj.mndoDerivedParameterRho(3) = 0.5/0.6685771472;
            obj.am1CoreintegralS = -52.028658 * Parameters.GetInstance().GetEV2AU();
            obj.am1CoreintegralP = -39.614239 * Parameters.GetInstance().GetEV2AU();
            obj.am1OrbitalExponentS = 1.808665;
            obj.am1OrbitalExponentP = 1.685116;
            obj.am1BondingParameterS = -15.715783 * Parameters.GetInstance().GetEV2AU();
            obj.am1BondingParameterP = -7.719283  * Parameters.GetInstance().GetEV2AU();
            obj.am1Alpha = 2.648274 / Parameters.GetInstance().GetAngstrom2AU();
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
            obj.am1ParameterK(1) = 0.011355 * Parameters.GetInstance().GetEV2AU();
            obj.am1ParameterK(2) = 0.045924 * Parameters.GetInstance().GetEV2AU();
            obj.am1ParameterK(3) =-0.020061 * Parameters.GetInstance().GetEV2AU();
            obj.am1ParameterK(4) =-0.001260 * Parameters.GetInstance().GetEV2AU();
            obj.am1ParameterL(1) = 5.00 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(2) = 5.00 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(3) = 5.00 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterL(4) = 5.00 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.am1ParameterM(1) = 1.60 * Parameters.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(2) = 1.85 * Parameters.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(3) = 2.05 * Parameters.GetInstance().GetAngstrom2AU();
            obj.am1ParameterM(4) = 2.65 * Parameters.GetInstance().GetAngstrom2AU();
            obj.am1DCoreintegralS = -52.183798 * Parameters.GetInstance().GetEV2AU();
            obj.am1DCoreintegralP = -39.368413 * Parameters.GetInstance().GetEV2AU();
            obj.am1DBondingParameterS = -15.682341 * Parameters.GetInstance().GetEV2AU();
            obj.am1DBondingParameterP = -7.804762 * Parameters.GetInstance().GetEV2AU();
            obj.am1DAlpha = 2.625506 / Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3CoreintegralS = -47.270320 * Parameters.GetInstance().GetEV2AU();
            obj.pm3CoreintegralP = -36.266918 * Parameters.GetInstance().GetEV2AU();
            obj.pm3OrbitalExponentS = 1.565085;
            obj.pm3OrbitalExponentP = 1.842345;
            obj.pm3BondingParameterS = -11.910015 * Parameters.GetInstance().GetEV2AU();
            obj.pm3BondingParameterP = -9.802755 * Parameters.GetInstance().GetEV2AU();
            obj.pm3Alpha = 2.707807 / Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3DerivedParameterD(1) = 0.0;
            obj.pm3DerivedParameterD(2) = 0.8332396384;
            obj.pm3DerivedParameterD(3) = 0.6647749859;
            obj.pm3DerivedParameterRho(1) = 0.5/0.4116151543;
            obj.pm3DerivedParameterRho(2) = 0.5/0.5885706542;
            obj.pm3DerivedParameterRho(3) = 0.5/0.7647513703;
            obj.pm3ParameterK(1) = 0.050107 * Parameters.GetInstance().GetEV2AU();
            obj.pm3ParameterK(2) = 0.050733 * Parameters.GetInstance().GetEV2AU();
            obj.pm3ParameterK(3) = 0.0;
            obj.pm3ParameterK(4) = 0.0;
            obj.pm3ParameterL(1) = 6.003165 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(2) = 6.002979 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3ParameterL(3) = 0.00;
            obj.pm3ParameterL(4) = 0.00;
            obj.pm3ParameterM(1) = 1.642214 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(2) = 0.892488 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3ParameterM(3) = 0.00;
            obj.pm3ParameterM(4) = 0.00;
            obj.pm3Gss = 11.200708 * Parameters.GetInstance().GetEV2AU();
            obj.pm3Gpp = 10.796292 * Parameters.GetInstance().GetEV2AU();
            obj.pm3Gsp = 10.265027 * Parameters.GetInstance().GetEV2AU();
            obj.pm3Gpp2 = 9.042566 * Parameters.GetInstance().GetEV2AU();
            obj.pm3Hsp = 2.290980 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralS = -48.241241 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgCoreintegralP = -36.461256 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgOrbitalExponentS = 1.567864;
            obj.pm3PddgOrbitalExponentP = 1.846659;
            obj.pm3PddgBondingParameterS = -11.952818 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgBondingParameterP =  -9.922411 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgAlpha = 2.725772 / Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3PddgDerivedParameterD(1) = 0.0;
            obj.pm3PddgDerivedParameterD(2) = 0.831413;
            obj.pm3PddgDerivedParameterD(3) = 0.663222;
            obj.pm3PddgDerivedParameterRho(1) = 1.214657;
            obj.pm3PddgDerivedParameterRho(2) = 0.848467;
            obj.pm3PddgDerivedParameterRho(3) = 0.652785;
            obj.pm3PddgParameterK(1) = 0.048906 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(2) = 0.047697 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgParameterK(3) = 0.0;
            obj.pm3PddgParameterK(4) = 0.0;
            obj.pm3PddgParameterL(1) = 5.765340 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(2) = 5.973721 / power(Parameters.GetInstance().GetAngstrom2AU(),2.0);
            obj.pm3PddgParameterL(3) = 0.00;
            obj.pm3PddgParameterL(4) = 0.00;
            obj.pm3PddgParameterM(1) = 1.682232 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(2) = 0.894406 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterM(3) = 0.00;
            obj.pm3PddgParameterM(4) = 0.00;
            obj.pm3PddgParameterPa(1) =-0.000743 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgParameterPa(2) = 0.000985 * Parameters.GetInstance().GetEV2AU();
            obj.pm3PddgParameterDa(1) = 0.836915 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3PddgParameterDa(2) = 1.585236 * Parameters.GetInstance().GetAngstrom2AU();
            obj.pm3DCoreintegralS = -47.275431 * Parameters.GetInstance().GetEV2AU();
            obj.pm3DCoreintegralP = -36.268916 * Parameters.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterS = -11.941466 * Parameters.GetInstance().GetEV2AU();
            obj.pm3DBondingParameterP = -9.819760 * Parameters.GetInstance().GetEV2AU();
            obj.pm3DAlpha = 2.721152 / Parameters.GetInstance().GetAngstrom2AU();
        end
        
    end
    
end