// Read field bounds
dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

dimensionedScalar rhoMin("rhoMin", rho.dimensions(), 0);
dimensionedScalar rhoMax("rhoMax", rho.dimensions(), 0);
fieldBounds.lookup("rho") >> rhoMin.value() >> rhoMax.value() ;

dimensionedScalar eMin("eMin", e.dimensions(), 0);
dimensionedScalar eMax("eMax", e.dimensions(), 0);
fieldBounds.lookup("e") >> eMin.value() >> eMax.value() ;


dimensionedScalar rhoEMin("rhoEMin", rhoE.dimensions(), 0);
dimensionedScalar rhoEMax("rhoEMax", rhoE.dimensions(), 0);
fieldBounds.lookup("rhoE") >> rhoEMin.value() >> rhoEMax.value();

dimensionedScalar UMax("UMax", U.dimensions(), 0);
fieldBounds.lookup("U") >> UMax.value();
