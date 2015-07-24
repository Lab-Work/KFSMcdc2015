package triggerer;

import org.jblas.DoubleMatrix;
import filters.Filter;

//Deterministic threshold-based sensor scheduler
public class DT extends Triggerer {	

	
	public DT() {
	}
	
	
	public void getNewParameters() {
	}	
	
	public void setNewParameters(DoubleMatrix _Z, double _zeta) {
        Z=_Z;
        zeta=_zeta;	
	}
	
	public void transmission(Filter _filter){
	    DoubleMatrix innov=_filter.measurements.sub(_filter.measure.mmul(_filter.f_mean));
	    double G=Math.sqrt(((_filter.measure.mmul(_filter.f_var).mmul(_filter.measure.transpose())).add(_filter.measureVar)).get(0,0));
	    phi=(1d / G)*(innov.get(0,0));
	    if (phi<=-zeta||phi>=zeta){
	    	gamma=1;
	    }
	    else{
	    	gamma=0;
	    }
	}
		
	
}