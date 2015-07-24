package filters;

import org.jblas.DoubleMatrix;
import doubleMatrix.GaussianGenerator;
import doubleMatrix.InverseMatrix;


public class KFSM extends Filter{
	
	
	public KFSM () {

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
			GaussianGenerator measureGenerator= new GaussianGenerator(DoubleMatrix.eye(sizeMeasurements).mul(0.25));
			DoubleMatrix noiseGenerator = measureGenerator.sample();	
			while (noiseGenerator.get(0,0)<=-triggerer.zeta || noiseGenerator.get(0,0)>triggerer.zeta){
				noiseGenerator = measureGenerator.sample();
			}
			double G=Math.sqrt(((measure.mmul(f_var).mmul(measure.transpose())).add(measureVar)).get(0,0));
			double innov=G*(noiseGenerator.get(0,0));
			DoubleMatrix innovation=DoubleMatrix.eye(sizeMeasurements);
			innovation.put(0,0,innov);
			
				
			Kgain = f_var.mmul(measure.transpose()).mmul(A).mul(beta);
			
			
            mean=f_mean.add(Kgain.mmul(innovation));	
            var=f_var.sub(Kgain.mmul(measure.mmul(f_var)));
		}
	}
		
	
}
