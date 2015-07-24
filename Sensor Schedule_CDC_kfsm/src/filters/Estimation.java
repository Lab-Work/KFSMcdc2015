package filters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import triggerer.*;
import model.*;
import org.jblas.DoubleMatrix;
import trueSolution.TrueSolution;
import trueSolution.TrueinGround;
import targetGroundTruth.TargetGroundTruth;

public class Estimation {
	

	Filter filter;

	public Estimation(String _nameFilter, Triggerer _triggerer, ModelKF _modelKF) {
		Filter f = Filter.createFilter(_nameFilter, _triggerer, _modelKF);	
		filter=f;
	}
		
	public void updateTrueSolution() {
		filter.modelKF.updateTrueSolution();
	}	

	public static void exportResult(TargetGroundTruth _targetGroundTruth, int limite, String folder) {
		try {
		    int numKFSM=100;
		    int numKFDT=100;
			int numFilters = numKFSM+numKFDT;
			int cells=_targetGroundTruth.cells;	
			
		    TrueSolution[] trueSolution=new TrueSolution[numFilters];
		    for(int i=0;i<numFilters;i++){
		        trueSolution[i] = new TrueinGround(_targetGroundTruth,cells,i, numFilters);
		    }
		    
			ModelKF [] modelKFs=new ModelKF[numFilters];
			for (int i=0;i<numFilters;i++){					
				modelKFs[i]=trueSolution[i].setModels();
			}
					
			BufferedWriter writerTrue;
			BufferedWriter writerRate;
			BufferedWriter writerTrVar;
			BufferedWriter[] writerFilter = new BufferedWriter[2];
			BufferedWriter[] writerError= new BufferedWriter[2];
			
			Triggerer []DTtriggerSM=new Triggerer [numKFSM];
			Triggerer []DTtriggerDT=new Triggerer [numKFDT];
			for(int i=0; i<numKFSM;i++){
				DTtriggerSM[i]=Triggerer.createTriggerer("DT");
				DTtriggerSM[i].setNewParameters(DoubleMatrix.eye(1).mul(200000), 0.70);
			}
			for(int i=0; i<numKFSM;i++){
				DTtriggerDT[i]=Triggerer.createTriggerer("DT");
				DTtriggerDT[i].setNewParameters(DoubleMatrix.eye(1).mul(200000), 0.70);
			}
						
			System.out.print("Beginning of simulation");
				
				new File("results/"+folder).mkdir();
				
				Estimation[] estimations = new Estimation[numFilters];
				double [] rate =new double [2];
				
				writerTrue = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"trueState.csv")));
				writerRate = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Rate.csv")));
				writerTrVar = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Trace.csv")));
				
				for (int i = 0; i<numFilters; i++) {
					if (i<numKFSM){
						estimations[i] = new Estimation("KFSM", DTtriggerSM[i], modelKFs[i]);
					}
					else {
						estimations[i] = new Estimation("KFDT", DTtriggerDT[i-numKFSM], modelKFs[i]);
					}

				}

				for (int i = 0; i<2; i++) {				
					String S = "results/"+folder+"/"+estimations[numKFSM-1+i].filter.getClass().getSimpleName();
					writerFilter[i] = new BufferedWriter(new FileWriter(new File(S+".csv")));
					writerError[i]=new BufferedWriter(new FileWriter(new File(S+"error.csv")));		
				}     
				
				for (int k=0; k<limite; k++) {
					
					for (int i=0;i<_targetGroundTruth.cells;i++){
						writerTrue.write(_targetGroundTruth.trueStatesGround.get(i,0)+",");
					}
					writerTrue.write("\n");		

					DoubleMatrix[] mean = new DoubleMatrix[2];	

					for (int i = 0; i<2; i++) {
						mean[i]= DoubleMatrix.zeros(_targetGroundTruth.trueStatesGround.getRows(), 1);
					}

					for (int i=0; i<2; i++) {
						if (i==0){
							for (int j=0;j<numKFSM; j++){
								mean[i]=mean[i].add((estimations[j].filter.mean).mul(((double)(1)/(double)(numKFSM))));
							}
						}
						else{
							for (int j=0;j<numKFDT; j++){
								mean[i]=mean[i].add((estimations[j+numKFSM].filter.mean).mul(((double)(1)/(double)(numKFDT))));
							}
						}
																					
						for (int j=0; j<mean[i].length ; j++) {
							writerFilter[i].write(mean[i].get(j)+",");							
						}
						
    					for (int j=0; j<mean[i].length ; j++) {
							writerError[i].write(mean[i].get(j)-_targetGroundTruth.trueStatesGround.get(j)+",");						
						}
    					writerError[i].write(Math.sqrt((((mean[i].sub(_targetGroundTruth.trueStatesGround)).transpose()).mmul(mean[i].sub(_targetGroundTruth.trueStatesGround))).get(0,0))+",");

						writerError[i].write("\n");
						writerFilter[i].write("\n");
					}
					writerTrVar.write("\n");
					
					double [] ratesep =new double [numFilters];
					for (int i=0; i<numFilters; i++) {
						ratesep[i]=ratesep[i]+((double)(estimations[i].filter.triggerer.gamma)/(double)(limite));
					}
	//				System.out.print(estimations[numKFSM].filter.triggerer.gamma);
					for (int i=0; i<2; i++) {
						if (i==0){
							for (int j=0;j<numKFSM; j++){
								rate[i]=rate[i]+((((double)(ratesep[j])/(double)(numKFSM))));
							}
						}
						else{
							for (int j=0;j<numKFDT; j++){
								rate[i]=rate[i]+((((double)(ratesep[j+numKFSM])/(double)(numKFDT))));
							}							
						}
					}
													
					for (int i=0; i<numFilters; i++) {
						estimations[i].filter.getNewParametersFromModel();
					}			
				
					_targetGroundTruth.update();			

					for(int i=0; i<numFilters; i++){
						trueSolution[i].update();
					}
					
					for(int i=0; i<numFilters; i++){
						trueSolution[i].updatemeasurement();
					}

					for (int i=0; i<numFilters; i++) {
						estimations[i].filter.NextStep();
					}					
		
					if (k==limite-1){
						writerRate.write(rate[0]+",");
						writerRate.write(rate[1]+",");
					}
				}
				
				for (int i = 0; i<2; i++) {
					writerFilter[i].flush(); writerFilter[i].close();
					writerError[i].flush(); writerError[i].close();
				}

				writerTrue.flush(); writerTrue.close();
				writerRate.flush(); writerRate.close();
				writerTrVar.flush(); writerTrVar.close();
				System.out.print(rate[0]);
				System.out.print(rate[1]);
			System.out.println(" End");			
		}
		catch (Exception e) {e.printStackTrace();}
		
	}
	
}
