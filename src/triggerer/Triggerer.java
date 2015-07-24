package triggerer;

import org.jblas.DoubleMatrix;
import filters.Filter;

public abstract class Triggerer {


	public int gamma;
	public DoubleMatrix Z;
	public double zeta;
	public double phi;
	
	abstract public void getNewParameters();
	abstract public void setNewParameters(DoubleMatrix _Z, double _DT);
	abstract public void transmission(Filter _filter);
			
 	//!!!!!!
 	
    public static Triggerer createTriggerer(String nameTriggerer) {
		String name = "triggerer."+nameTriggerer;
 		try {
			@SuppressWarnings("rawtypes")
			Class cl = Class.forName(name);
			@SuppressWarnings("unchecked")
			Triggerer trig = (Triggerer) cl.getConstructor().newInstance();		
			trig.Z=DoubleMatrix.eye(1).mul(40000);
			trig.zeta=0.006;
			return trig;
 		}
		catch(Throwable t) {System.out.println(t);}
				
		return null;
	}
	

}
	

