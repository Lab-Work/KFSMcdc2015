package trueSolution;
import targetGroundTruth.TargetGroundTruth;;

public class TrueinGround extends TrueSolution{
	
	public TrueinGround (TargetGroundTruth _trueGround,int _cells,  int _index, int _numSecs){
		initial(_trueGround,_cells,_index,_numSecs);
	}

	public void propagate() {
	trueStates=trueGround.trueStatesGround.dup();	
	}
	
}
