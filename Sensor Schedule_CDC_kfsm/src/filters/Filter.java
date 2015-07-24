package filters;

import model.ModelKF;
import org.jblas.DoubleMatrix;
import triggerer.*;

public abstract class Filter {

	public DoubleMatrix var;
	public DoubleMatrix varUpper;
	public DoubleMatrix varLower;
	public DoubleMatrix mean;
	public DoubleMatrix f_var;
	public DoubleMatrix f_varUpper;
	public DoubleMatrix f_varLower;
	public DoubleMatrix f_mean;		
	public DoubleMatrix Kgain;	
	public DoubleMatrix measure;
	public DoubleMatrix measurements;
	public DoubleMatrix measureVar;	
	public double rate;
	public double beta;
	
	
	int size;
	int sizeMeasurements;
	int cells;	
	int stepMeasurements;
	int numUp;
	
	DoubleMatrix modelVar;
	DoubleMatrix priorVar;
	
	Triggerer triggerer;
	ModelKF modelKF;
	
	
	
	abstract public void NextStep();

		
	public void initial(Triggerer _triggerer, ModelKF _modelKF) {
		modelKF=_modelKF;
		triggerer=_triggerer;	
		getNewParametersFromModel();		
		mean = modelKF.initialMean;
		f_mean=mean.dup();
		var =modelKF.initialVar;
		f_var=var.dup();
		varUpper=var.dup();
		varLower=var.dup();
		f_varUpper=var.dup();
		f_varLower=var.dup();
		priorVar=var.dup();		
		numUp = 0;
		double betaU=2*(triggerer.zeta)*(Math.exp(-(triggerer.zeta)*(triggerer.zeta) / 2) / Math.sqrt(2 * Math.PI));
		double betaL=1-(2*Q(triggerer.zeta));
		beta=betaU / betaL;
		rate=1-betaL;
	}
	

    
	public static Filter createFilter(String name,Triggerer _triggerer, ModelKF _modelKF) {
		Filter f = null;
		if (name.equals("KFSM")) f = new KFSM();
		else if (name.equals("KFDT")) f = new KFDT();
		else if (name.equals("KFST")) f = new KFST();
		else throw new Error("Filter not recognized");
		f.initial(_triggerer, _modelKF);
		return f;	 	
	}
	

	
	protected DoubleMatrix getMeasurements() {
		return modelKF.getMeasureVector();
	}
	
	protected void computeTriggerer() {
		triggerer.transmission(this);
	}
	
	protected DoubleMatrix computeVar(DoubleMatrix _Var){
		DoubleMatrix _fVar=DoubleMatrix.zeros(_Var.getRows(), _Var.getColumns());
		_fVar=modelKF.trueSolution.AMatrix.mmul(_Var).mmul(modelKF.trueSolution.AMatrix.transpose()).add(modelVar);
		return _fVar;
	}
	
	protected void getNewParametersFromModel() {
		size = modelKF.size;
		modelVar = modelKF.modelVarKF;
		sizeMeasurements = modelKF.sizeMeasurements;
		stepMeasurements = modelKF.stepMeasurements;
		measureVar = modelKF.measureVar;
		measure = modelKF.measure;
	}
	
    protected static double phi(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }
	
	protected static double Q (double z){
		if (z < -8.0) return 0.0;
        if (z >  8.0) return 1.0;
        double sum = 0.0, term = z;
        for (int i = 3; sum + term != sum; i += 2) {
            sum  = sum + term;
            term = term * z * z / i;
        }
        return 1-(0.5 + sum * phi(z));
        
	}
	
}
