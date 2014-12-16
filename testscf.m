function testscf()


import SEQC.*;

%%
disp('Begin cndo/2 tests');

c1 = AtomFactory.Create(EnumAtom.C, 1, [         0         0         0]./0.5291772); % ang 2 au
h2 = AtomFactory.Create(EnumAtom.Cl, 2, [         0         0    1.6822]./0.5291772);
h3 = AtomFactory.Create(EnumAtom.H, 3, [    1.0204         0   -0.3602]./0.5291772);
h4 = AtomFactory.Create(EnumAtom.H, 4, [   -0.5102   -0.8835   -0.3604]./0.5291772);
h5 = AtomFactory.Create(EnumAtom.H, 5, [   -0.5100    0.8836   -0.3604]./0.5291772);

mol = Molecule();
mol.AddAtom(c1);
mol.AddAtom(h2);
mol.AddAtom(h3);
mol.AddAtom(h4);
mol.AddAtom(h5);
mol.CalcBasics();

model = Cndo2();
model.SetMolecule(mol);

tic
model.DoSCF();
toc

epsRel = 1e-4;
conditionEScf = testRelErr(model.elecSCFEnergy, -2.552489e+01) < epsRel;
conditionEMo = (testRelErr(model.energiesMO(7), -5.142121e-01) < epsRel)...
    && (testRelErr(model.energiesMO(8), 1.173577e-01) < epsRel);
conditionCoreRep = testRelErr(model.coreRepulsionEnergy, 2.044227e+01) < epsRel;
dip = reshape(model.electronicTransitionDipoleMoments, 3, []);
conditionDipole = (testRelErr(dip(1), -3.114910e-04) < epsRel)...
    && (testRelErr(dip(2), -2.750244e-04) < epsRel)...
    && (testRelErr(dip(3), 2.206855e+01) < epsRel);
if(conditionEScf && conditionEMo && conditionCoreRep && conditionDipole)
    disp('cndo/2 tests all passed');
else
    disp('some cndo/2 tests failed');
end
if(~conditionEScf)
    disp('cndo/2 scf energy failed');
end
if(~conditionEMo)
    disp('cndo/2 homo or lumo failed');
end
if(~conditionCoreRep)
    disp('cndo/2 core repulsion energy failed')
end
if(~conditionDipole)
    disp('cndo/2 dipole failed')
end

%%
disp('Begin indo tests');

c1 = AtomFactory.Create(EnumAtom.C, 1, [         0         0         0]./0.5291772); % ang 2 au
h2 = AtomFactory.Create(EnumAtom.F, 2, [         0         0    1.6822]./0.5291772);
h3 = AtomFactory.Create(EnumAtom.H, 3, [    1.0204         0   -0.3602]./0.5291772);
h4 = AtomFactory.Create(EnumAtom.H, 4, [   -0.5102   -0.8835   -0.3604]./0.5291772);
h5 = AtomFactory.Create(EnumAtom.H, 5, [   -0.5100    0.8836   -0.3604]./0.5291772);

mol = Molecule();
mol.AddAtom(c1);
mol.AddAtom(h2);
mol.AddAtom(h3);
mol.AddAtom(h4);
mol.AddAtom(h5);
mol.CalcBasics();

model = Indo();
model.SetMolecule(mol);

tic
model.DoSCF();
toc

epsRel = 1e-4;
conditionEScf = testRelErr(model.elecSCFEnergy, -3.536833e+01) < epsRel;
conditionEMo = (testRelErr(model.energiesMO(7), -6.572898e-01) < epsRel)...
    && (testRelErr(model.energiesMO(8), 1.104153e-01) < epsRel);
conditionCoreRep = testRelErr(model.coreRepulsionEnergy, 2.044227e+01) < epsRel;
dip = reshape(model.electronicTransitionDipoleMoments, 3, []);
conditionDipole = (testRelErr(dip(1), -3.355399e-04) < epsRel)...
    && (testRelErr(dip(2), -2.566653e-04) < epsRel)...
    && (testRelErr(dip(3), 1.170728e+01) < epsRel);
if(conditionEScf && conditionEMo && conditionCoreRep && conditionDipole)
    disp('indo tests all passed');
else
    disp('some indo tests failed');
end
if(~conditionEScf)
    disp('indo scf energy failed');
end
if(~conditionEMo)
    disp('indo homo or lumo failed');
end
if(~conditionCoreRep)
    disp('indo core repulsion energy failed')
end
if(~conditionDipole)
    disp('indo dipole failed')
end



%%
disp('Begin mndo tests');
c1 = AtomFactory.Create(EnumAtom.C, 1, [         0         0         0]./0.5291772); % ang 2 au
h2 = AtomFactory.Create(EnumAtom.Cl, 2, [         0         0    1.6822]./0.5291772);
h3 = AtomFactory.Create(EnumAtom.H, 3, [    1.0204         0   -0.3602]./0.5291772);
h4 = AtomFactory.Create(EnumAtom.H, 4, [   -0.5102   -0.8835   -0.3604]./0.5291772);
h5 = AtomFactory.Create(EnumAtom.H, 5, [   -0.5100    0.8836   -0.3604]./0.5291772);

mol = Molecule();
mol.AddAtom(c1);
mol.AddAtom(h2);
mol.AddAtom(h3);
mol.AddAtom(h4);
mol.AddAtom(h5);
mol.CalcBasics();

model = Mndo();
model.SetMolecule(mol);

tic
model.DoSCF();
toc

epsRel = 1e-4;
conditionEScf = testRelErr(model.elecSCFEnergy, -1.931167e+01) < epsRel;
conditionEMo = (testRelErr(model.energiesMO(7), -4.481019e-01) < epsRel)...
    && (testRelErr(model.energiesMO(8), 5.382423e-02) < epsRel);
conditionCoreRep = testRelErr(model.coreRepulsionEnergy, 1.740992e+01) < epsRel;
dip = reshape(model.electronicTransitionDipoleMoments, 3, []);
conditionDipole = (testRelErr(dip(1), -3.189724e-04) < epsRel)...
    && (testRelErr(dip(2), -2.759633e-04) < epsRel)...
    && (testRelErr(dip(3), 2.027049e+01) < epsRel);
if(conditionEScf && conditionEMo && conditionCoreRep && conditionDipole)
    disp('mndo tests all passed');
else
    disp('some mndo tests failed');
end
if(~conditionEScf)
    disp('mndo scf energy failed');
end
if(~conditionEMo)
    disp('mndo homo or lumo failed');
end
if(~conditionCoreRep)
    disp('mndo core repulsion energy failed')
end
if(~conditionDipole)
    disp('mndo dipole failed')
end

%%
disp('Begin am1 tests');

model = Am1();
model.SetMolecule(mol);

tic
model.DoSCF();
toc

epsRel = 1e-4;
conditionEScf = testRelErr(model.elecSCFEnergy, -1.996811e+01) < epsRel;
conditionEMo = (testRelErr(model.energiesMO(7), -4.153534e-01) < epsRel)...
    && (testRelErr(model.energiesMO(8), 6.685816e-02) < epsRel);
conditionCoreRep = testRelErr(model.coreRepulsionEnergy, 1.733616e+01) < epsRel;
dip = reshape(model.electronicTransitionDipoleMoments, 3, []);
conditionDipole = (testRelErr(dip(1), -3.054159e-04) < epsRel)...
    && (testRelErr(dip(2), -2.695722e-04) < epsRel)...
    && (testRelErr(dip(3), 2.089521e+01) < epsRel);
if(conditionEScf && conditionEMo && conditionCoreRep && conditionDipole)
    disp('am1 tests all passed');
else
    disp('some am1 tests failed');
end
if(~conditionEScf)
    disp('am1 scf energy failed');
end
if(~conditionEMo)
    disp('am1 homo or lumo failed');
end
if(~conditionCoreRep)
    disp('am1 core repulsion energy failed')
end
if(~conditionDipole)
    disp('am1 dipole failed')
end


%%
disp('Begin pm3 tests');

model = Pm3();
model.SetMolecule(mol);

tic
model.DoSCF();
toc

epsRel = 1e-4;
conditionEScf = testRelErr(model.elecSCFEnergy, -1.769142e+01) < epsRel;
conditionEMo = (testRelErr(model.energiesMO(7), -3.851521e-01) < epsRel)...
    && (testRelErr(model.energiesMO(8), 7.160589e-02) < epsRel);
conditionCoreRep = testRelErr(model.coreRepulsionEnergy, 1.736127e+01) < epsRel;
dip = reshape(model.electronicTransitionDipoleMoments, 3, []);
conditionDipole = (testRelErr(dip(1), -3.188360e-04) < epsRel)...
    && (testRelErr(dip(2), -2.665222e-04) < epsRel)...
    && (testRelErr(dip(3), 2.078154e+01) < epsRel);
if(conditionEScf && conditionEMo && conditionCoreRep && conditionDipole)
    disp('pm3 tests all passed');
else
    disp('some pm3 tests failed');
end
if(~conditionEScf)
    disp('pm3 scf energy failed');
end
if(~conditionEMo)
    disp('pm3 homo or lumo failed');
end
if(~conditionCoreRep)
    disp('pm3 core repulsion energy failed')
end
if(~conditionDipole)
    disp('pm3 dipole failed')
end

end

function res = testRelErr(numTar, numRef)
res = abs((numTar - numRef) / numRef);
end
