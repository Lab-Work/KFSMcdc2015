package trueSolution;

import org.jblas.DoubleMatrix;
import doubleMatrix.GaussianGenerator;
import targetGroundTruth.TargetGroundTruth;
import model.*;



public abstract class TrueSolution {

	public TargetGroundTruth trueGround;
	
	public int cells;
	public int index;
	
	public int numSections;
		
	public DoubleMatrix measureVarDensity;//the measurement error covariance matrix
	public GaussianGenerator measureGeneratorDensity;	
	public DoubleMatrix measurementsDensity;//sensor data
	
	public DoubleMatrix trueStates;
	public DoubleMatrix currenttrueStates;	//used for computing evolution of true density
	
	public int stepMeasurements=1;//number of steps between two measurements
	public int sizeMeasurements; //for one state vector

	public int numUp;
	
	public DoubleMatrix AMatrix;
	

	abstract public void propagate();

	
	public void initial(TargetGroundTruth _trueGround, int _cells, int _index, int _numSecs) {
		
		trueGround=_trueGround;
		
		index=_index;
		
		cells=_cells;
		numSections=_numSecs;	
				
		sizeMeasurements=1;//**upstream and downstream measurements
		//**measurement error covariance matrix, which is larger than the true measurement error 
		//to ensure satisfactory performance of the filter
		measureVarDensity = DoubleMatrix.eye(sizeMeasurements).mul(1);
		measureGeneratorDensity = new GaussianGenerator(measureVarDensity.mul(0.01));

		trueStates=DoubleMatrix.zeros(cells,1);
		currenttrueStates=trueStates.dup();

		//**Start setting initial condition for true solution		
		
		trueStates=trueGround.trueStatesGround.dup();
		//**end setting initial condition of the true state
		
		AMatrix=trueGround.AMatrix;
		numUp = 0;
	}


	public void update() {		
		numUp++;
		propagate();		
	}
	
    public void updatemeasurement() {	
		newMeasurements();		
	}
	
	private void newMeasurements() {
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) getMeasurementsAll();
	}
	
	private void getMeasurementsAll() {
		
		measurementsDensity = DoubleMatrix.zeros(sizeMeasurements, 1);
		DoubleMatrix noiseDensity = measureGeneratorDensity.sample();		
		for (int i=0; i<sizeMeasurements; i++) {			
			double pointDensity = trueStates.get((i+1),0);
//			System.out.print(i+1);
//			double pointDensity = trueStates.get(i*distMeasurements,0);
			pointDensity+=noiseDensity.get(i);//comment this line to turn off measurement noise
			measurementsDensity.put(i, pointDensity);
		}		
	}
	
 	public ModelKF setModels(){
 		ModelKF ms;
 		ms=ModelKF.createModel(this,"filterModel");		
 		return ms;
 	}


}
