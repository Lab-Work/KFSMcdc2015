package filters;

import org.jblas.DoubleMatrix;
import doubleMatrix.InverseMatrix;

public class KFST extends Filter{
	
	
	public KFST () {

	}
	
	public void NextStep() {
		//getNewParametersFromModel();
		//measureGenerator = new GaussianGenerator(measureVar);		
		forecast();
		analysis();
		numUp++;
	}
	
	private void forecast() {
		f_mean = modelKF.propagate(mean);
		f_var = computeVar(var);
		f_varUpper=computeVar(varUpper);
		f_varLower=computeVar(varLower);
		priorVar=var;		
	}

	
	private void analysis() {
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) {
			measurements = getMeasurements();
			
			computeTriggerer();
			
			updateState();
		}
		else {mean =f_mean.dup(); var = f_var.dup();}
//		mean = computeMean(samples); var = f_var;
	}
	
	
	
	
	
	private void updateState() {
		DoubleMatrix A=new DoubleMatrix(modelKF.sizeMeasurements,modelKF.sizeMeasurements);
		
		if (triggerer.gamma==1){
			A = InverseMatrix.invPoSym( (measure.mmul(f_var.mmul(measure.transpose()))).add(measureVar) );
			Kgain = f_var.mmul(measure.transpose()).mmul(A);		
//			measureVector.print();
			mean=f_mean.add(Kgain.mmul(measurements.sub(measure.mmul(f_mean))));	
			var = f_var.sub(Kgain.mmul(measure.mmul(f_var)));
		}
		else{
			A = InverseMatrix.invPoSym((measure.mmul(f_var.mmul(measure.transpose()))).add(measureVar).add(InverseMatrix.invPoSym(triggerer.Z)));
			
			Kgain = f_var.mmul(measure.transpose()).mmul(A);
            mean=f_mean.dup();
            var=f_var.sub(Kgain.mmul(measure.mmul(f_var)));
		}
		
		DoubleMatrix AU=new DoubleMatrix(modelKF.sizeMeasurements,modelKF.sizeMeasurements);
		DoubleMatrix AL=new DoubleMatrix(modelKF.sizeMeasurements,modelKF.sizeMeasurements);
		AU = InverseMatrix.invPoSym( (measure.mmul(f_varUpper.mmul(measure.transpose()))).add(measureVar).add(InverseMatrix.invPoSym(triggerer.Z)) );
		AL = InverseMatrix.invPoSym( (measure.mmul(f_varLower.mmul(measure.transpose()))).add(measureVar) );
		
		varUpper = f_varUpper.sub((f_varUpper.mmul(measure.transpose()).mmul(AU)).mmul(measure.mmul(f_varUpper)));
		varLower = f_varLower.sub((f_varLower.mmul(measure.transpose()).mmul(AL)).mmul(measure.mmul(f_varLower)));
		
	}
	

		
	
}
