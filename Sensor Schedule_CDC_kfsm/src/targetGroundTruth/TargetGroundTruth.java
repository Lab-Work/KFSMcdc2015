package targetGroundTruth;

import org.jblas.DoubleMatrix;

public abstract class TargetGroundTruth {

	public int cells;
	
	public DoubleMatrix trueStatesGround;
	public DoubleMatrix trueStatesGroundPrior;
	
	public int numUp;
	
	public DoubleMatrix AMatrix;
		
	abstract public void propagateGround();
	
	public void initial(int _cells) {
	
		cells=_cells;

		trueStatesGround=DoubleMatrix.zeros(cells,1);
		
		//**Start setting initial condition for the true state
	    trueStatesGround.put(0,0,2);
	    trueStatesGround.put(1,0,1);
		//**End setting initial condition for the true state
	    trueStatesGroundPrior=trueStatesGround.dup();
	    
	    //**Start setting A Matrix
	    AMatrix=DoubleMatrix.zeros(cells,cells);
	    AMatrix.put(0,0,1.05);
	    AMatrix.put(1,0,0.1);
	    AMatrix.put(1,1,0.9);
	    //**End setting A Matrix
	    		
		numUp = 0;
	}



	public void update() {		
		numUp++;
		propagateGround();	
	}
	

}
