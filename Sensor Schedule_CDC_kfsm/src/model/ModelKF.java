package model;

import org.jblas.DoubleMatrix;
import org.jfree.data.xy.XYSeries;
import doubleMatrix.Concat;
import doubleMatrix.GaussianGenerator;
import trueSolution.TrueinGround;
import trueSolution.TrueSolution;


@SuppressWarnings("unused")
public abstract class ModelKF {

	abstract public void initialThis();
	abstract public DoubleMatrix propagate(DoubleMatrix _density);
	abstract public DoubleMatrix getMeasureVector();	

	public DoubleMatrix initialMean;//initial estimate
	public DoubleMatrix initialVar;//initial estimation error covariance
	
	public int sizeMeasurements; 
	public int stepMeasurements;
	public int size; 
	
	public DoubleMatrix measureVar;//measurement error covariance matrix
	public DoubleMatrix measure;//output matrix

	public TrueSolution trueSolution;
	
	public String nameModel;
	public GaussianGenerator initialGeneratorDensity;
	
	int cells;//number of cells in a section
	
	DoubleMatrix initialGuess;	
	public DoubleMatrix modelVarKF;//used to set modelVar&modelVar in Filter
	DoubleMatrix initialGuessVar;//used to generate noise of initial estimate
	
	public void initial(TrueSolution _trueSolution) {
		
		trueSolution = _trueSolution;
		cells = trueSolution.cells;	
		stepMeasurements = trueSolution.stepMeasurements;

		defineVarAndMean();
		initialThis(); 
		
 		nameModel = getClass().getSimpleName()+"_";
		nameModel+=Integer.toString((int) Math.floor(trueSolution.index));
		nameModel+="_";
	}
		
	
	public static ModelKF createModel(TrueSolution _trueSolution, String nameEvolution ) {
 		String name = "model."+nameEvolution;
 		try {			
			@SuppressWarnings("rawtypes")
			Class cl = Class.forName(name);			
			@SuppressWarnings("unchecked")
			ModelKF m = (ModelKF) cl.getConstructor().newInstance(); 
			m.initial(_trueSolution);
			return m;
		}
		catch(Throwable t) {System.out.println("Model creation has failed");t.printStackTrace();}
		return null;
	}
 	
 	public void defineVarAndMean() {
 		
 		//**start setting initial mean and var
		DoubleMatrix initialGuess=DoubleMatrix.zeros(cells, 1);
		initialGuess.put(0,0,10);
		initialGuess.put(1,0,1.05);		
		
		initialMean=initialGuess;
		initialVar=DoubleMatrix.eye(cells).mul(100);
		//**end setting initial mean and var
		
		modelVarKF = DoubleMatrix.eye(cells).mul(0.01);//**model error covariance matrix
		modelVarKF.put(1,1,1);
 	}
 	
	
	public void updateTrueSolution() {
		trueSolution.update();
	}
	
	
	//*here update model error covariance matrix if needed
	public void updateModelVarKF(DoubleMatrix _ModelVar) {
        modelVarKF=_ModelVar;
	}	
}
