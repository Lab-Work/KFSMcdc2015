package targetGroundTruth;

import org.jblas.DoubleMatrix;

public class GroundTruth extends TargetGroundTruth{
	
	public GroundTruth(int _cells){
		initial(_cells);
	}

	public void propagateGround() {
		DoubleMatrix _densitynext=DoubleMatrix.zeros(trueStatesGround.getRows(), 1);
		trueStatesGroundPrior=trueStatesGround.dup();
		_densitynext=AMatrix.mmul(trueStatesGroundPrior);
		trueStatesGround=_densitynext.dup();
	}
	
}
