{
	#include "updateFusionFields.H"
	#include "updateVaporFields.H"
    radiation->correct();
    rhoCp =  rho*(fluid.Cp() + cp_fusion);

    const surfaceScalarField rhoCpPhi(fvc::interpolate(fluid.Cp() +  cp_fusion)*rhoPhi);

    const volScalarField kappaEff
    (
        "kappaEff",
        fluid.kappa() + fluid.Cp()*turbulence->mut()/fluid.Prt()
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(rhoCpPhi,T)")
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
      - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
      ==
        fluid.heatTransfer(T)
      + radiation->ST(T)
      + fvOptions(rhoCp, T)
	  + Q_vapor
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    fluid.correct();
	
	T = min(max(T, dimensionedScalar(dimTemperature, 298.0)), dimensionedScalar(dimTemperature, 5000.0));
	

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
		
		forAll(C, i)
		{
			if(molten[i] < 1.0)
			{
				if(beta[i] >= 0.0)
					molten[i] = beta[i];
			}
		}
		molten.correctBoundaryConditions();
}
