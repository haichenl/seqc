c1 = AtomFactory.Create(AtomType.C, 1, [         0         0         0]./0.5291772); % ang 2 au
h2 = AtomFactory.Create(AtomType.Cl, 2, [         0         0    1.6822]./0.5291772);
h3 = AtomFactory.Create(AtomType.H, 3, [    1.0204         0   -0.3602]./0.5291772);
h4 = AtomFactory.Create(AtomType.H, 4, [   -0.5102   -0.3835   -0.3604]./0.5291772);
h5 = AtomFactory.Create(AtomType.H, 5, [   -0.5100    0.3836   -0.3604]./0.5291772);

mol = Molecule();
mol.AddAtom(c1);
mol.AddAtom(h2);
mol.AddAtom(h3);
mol.AddAtom(h4);
mol.AddAtom(h5);
mol.CalcBasics();

cndo2 = Cndo2();
cndo2.SetMolecule(mol);

%%
cndo2.DoSCF();

