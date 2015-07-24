package filters;

import org.jblas.DoubleMatrix;
import doubleMatrix.InverseMatrix;


public class KFDT extends Filter{
	
	
	public KFDT () {
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
		A = InverseMatrix.invPoSym( (measure.mmul(f_var.mmul(measure.transpose()))).add(measureVar) );

		if (triggerer.gamma==1){
			Kgain = f_var.mmul(measure.transpose()).mmul(A);		
//			measureVector.print();
			mean=f_mean.add(Kgain.mmul(measurements.sub(measure.mmul(f_mean))));	
			var = f_var.sub(Kgain.mmul(measure.mmul(f_var)));
		}
		else{
			Kgain = f_var.mmul(measure.transpose()).mmul(A).mul(beta);
            mean=f_mean.dup();
            var=f_var.sub(Kgain.mmul(measure.mmul(f_var)));
		}
		
	}
	

		
	
}
