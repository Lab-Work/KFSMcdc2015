package triggerer;

import org.jblas.DoubleMatrix;
import filters.Filter;
import java.util.Random;

//Stochastic threshold-based sensor scheduler
public class ST extends Triggerer {	
	
	public ST() {
	}
	
	
	public void getNewParameters() {
	}	
	
	public void setNewParameters(DoubleMatrix _Z, double _zeta) {
        Z=_Z;
        zeta=_zeta;	
	}
	
	public void transmission(Filter _filter){
	    Random rnd = new Random();
	    double st=rnd.nextDouble();
	    DoubleMatrix innov=_filter.measurements.sub(_filter.measure.mmul(_filter.f_mean));
	    double phi=Math.exp(-0.5*(innov.transpose().mmul(Z).mmul(innov).get(0,0)));
	    if (st<=phi){
	    	gamma=0;
	    }
	    else{
	    	gamma=1;
	    }
	}
		
	
}