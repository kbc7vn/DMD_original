package muscle;


import java.util.ArrayList;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.query.space.grid.VNQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.SpatialMath;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;

/**
 * @author Kelley Virgilio
 *
 */

public class Fibroblast {
	
	// Fibroblast parameters
	private static Grid<Object> grid; 
	private static ContinuousSpace<Object> space;
	public static double eosinophilRecruit; // eosinphil recruitment parameter
	public static double fibroblastRecruit; // fibroblast recruitment
	public static double fibroblastApop; // fibroblast apoptosis 
	public static double fibroblastBlockApop; //fibroblast block apoptosis
	public static double fibroblastProlif; // fibroblast proliferation
	public static double fibroblastQuiescence; //fibroblast quiescence
	public static double fibroblastExpansion; //fibroblast expansion
	public static double fibrobRecruitSaturation; // fibroblast recruitment saturation limit
	public static double fibroExpSaturation; // fibroblast expansion saturation limit
	public static final int myofibroblastSwitch = 12; // time of switch from fibroblast to myofibroblast
	public static final int fibrobApopCount = 3; // fibroblast apoptosis time
	public static final int expTime = 12; // fibroblast expansion time
	public int phenotype; // fibroblast phenotype: starts at 0, 1+ converting to myofibroblast with more collagen
	private int recruitAge; // fibroblast recruitment time- counter
	private int expansionAge; // fibroblast expansion time- counter
	private int apopAge; // fibroblast apoptosis time- counter
	private int resident;// active fibroblast that is there at start of simulation
	
	// Test parameters:
	public static final double murphRep = 1; // 1 = control
	public static final double murphRepRecruitParam = 1; // 1 = control
	private static final int fibroblastSat = 10; // 10 control
	public static final double fibroScale = .5; // number of fibroblasts/fiber; .5 control
	
	// DISEASE STATE PARAMETERS: 
	// Murphy replication
//	public static final double murphRep = .4; // 1 = control; .4 at murphy sim
//	public static final double murphRepRecruitParam = .8; // 1 = control; .4 at murphy sim
//	public static final double fibroScale = .2; // number of fibroblasts/fiber
	// Fibroblast parameter analysis:
//	public static final double murphRepRecruitParam = 1.2; // 1 = control
//	public static final double fibroScale = 1.5; // number of fibroblasts/fiber
	//public static final double murphRep = .4; // 1 = control; .4 at murphy sim
	
	public static int count = 0;
	private int daughter; //same as ssc-- keeps track of the number of divisions
	
	public Fibroblast(Grid<Object> grid, ContinuousSpace<Object> space, int phenotype, int recruitAge, int apopAge, int resident, int expansionAge, int daughter)
	{
		this.grid = grid;
		this.space = space;
		this.phenotype = phenotype;// starts at 0, 1+ converting to myofibroblast with more collagen
		this.recruitAge = recruitAge; // fibroblast recruitment time- counter
		this.apopAge = apopAge; // fibroblast apoptosis time- counter
		this.resident = resident; // 0 = active fibroblast; 1 = resident/quiescent will restore pool of fibroblasts
		this.expansionAge = expansionAge; // fibroblast expansion time- counter
		this.daughter = daughter; // records number of divisions- similar to SSCs
	}
	
	//FIBROBLAST BEHAVIORS
	@ScheduledMethod(start = 2, interval = 1, priority = 2)
	public void fibroblastStep(){
		Context context = ContextUtils.getContext(this);
		// MOVE
		if (this.getPhenotype() < myofibroblastSwitch && this.resident == 0 && this.expansionAge == 0){
			move(); 
		}
		// MYOFIBROBLAST SECRETIONS
		if (this.phenotype >= myofibroblastSwitch){
			myofibroblastSecretions();
		}
		// ACTIVATION
		if (this.resident == 1 && this.getRecruitAge() == 0){ // 15 hour time period for recruitment
			// Check to see if the fibroblasts should be "activated"
			fibrobActivation();
		}
		// SENSE ENVIRONMENT
		if (this.resident == 0 && this.expansionAge == 0){
			senseEnvironment();
		}
		// RECRUITMENT -- called from Fiber
		// until the fibroblast has been fully recruited it is not counted as a fibroblast
		if (this.getRecruitAge() >= 1){
			setRecruitAge(this.getRecruitAge() + 1); // tracks fibroblast time to recruitment
		}
		if (this.getRecruitAge() >= 15 - RandomHelper.nextIntFromTo(0, 7)){ // 7-15 hours for recruitment
			this.resident = 0; // make it active
			this.setRecruitAge(0);
		}
		// EXPANSION/PROLIFERATION
		if (this.resident == 0 && this.getRecruitAge() == 0){
			fibrobExpansion();
		}
		if (this.expansionAge >= 1){
			this.setExpansionAge(this.expansionAge + 1);
		}
		// APOPTOSIS
		if (this.resident == 0 && this.expansionAge == 0 && this.getPhenotype() < myofibroblastSwitch){
			fibrobApoptosis();
		}
		if (this.getApopAge() >= fibrobApopCount && this.resident == 0){ // committed to apoptosis
			setApopAge(this.getApopAge() + 1);
		}
		// "QUIESCENCE" = HOMEOSTASIS/BASELINE CONDITIONS
		if (this.getPhenotype() < myofibroblastSwitch  && this.getRecruitAge() == 0 && this.getApopAge() < 2 && this.resident == 0 && this.expansionAge == 0){
			fibrobQuiescence();
		}
		
		// DISEASE STATE AND TEST PARAMETERS:
		// MURPHY REPLICATION-- blocking a % of fibroblasts
//		if (RandomHelper.nextIntFromTo(0, 50) < 1){
//			context.remove(this);
//		}
		// MACROPHAGE DEPLETION
		if (Fiber.macDepletion == 1 && InflamCell.getTick() > 100){
			// only used if trying to replicate macrophage depletion
			fibrobNecrosisRemove();
		}
	}
	
	// FIBROBLAST SENSE ENVIRONMENT
	public void senseEnvironment(){
		Context context = ContextUtils.getContext(this);
		double[] growthFactors = GrowthFactors.getGrowthFactors();
		// DISEASE STATE PARAMETERS
		double tnf = growthFactors[1] + Fiber.mdxTnf;
		double ifn = growthFactors[15] + Fiber.mdxIfn;
		// FIBROBLAST PROLIFERATION AND APOPTOSIS
		fibroblastRecruit = InflamCell.getFibrobRecruit(); //IL4-mediated eosinophil recruitment of NMMPs + SSC recruitment;
		fibroblastExpansion = InflamCell.inflamCells[8]; // SSC-dependent proliferation of fibroblasts (Murphy et al. 2011)
		fibroblastApop = tnf; // tnf induced fibroblast apotosis (Lemos 2015)
		fibroblastBlockApop = GrowthFactors.getActiveTgf(InflamCell.getTick()); //active tgf blocks fibroblast apotosis (Lemos 2015)
		fibroblastProlif = 0; // see fibroblastExpansion
		double fiberNorm = Fiber.origFiberNumber; // normalize number of fibroblasts to the number of fibers
		if (Fiber.origFiberNumber < getFibroblasts(context).size()){
			// If there are more fibroblasts than fibers-- then normalize by the number of fibroblasts
			fiberNorm = getFibroblasts(context).size();
		}
		// Saturation limits (most not utilized in typical simulations- mostly used at chronic parameter testing)
		fibrobRecruitSaturation = (fibroblastRecruit/fiberNorm)*15*Fiber.tcf4Scale; // saturation limits for fibroblast recruitment
		fibroExpSaturation = (1 - (Fiber.getFiberElems(context).size()/(Fiber.origCsaMean*Fiber.origFiberNumber)))*Fiber.origFiberNumber*10*Fiber.tcf4Scale; // Fibroblast expansion saturation
		fibroblastQuiescence = (Fiber.origCsaMean*Fiber.origFiberNumber - Fiber.getFiberElems(context).size())/(Fiber.origFiberNumber)*7; // returns fibroblasts to homeostatic condition
	}
	
	// FIBROBLAST ACTIVATION
	public void fibrobActivation(){
		// Activate resident fibroblasts if there is a lot of damage
		// All fibroblasts start off as 'resident'
		Context context = ContextUtils.getContext(this);
		if (this.resident == 1 && Necrosis.getPercentNecrotic(context) > .001 && RandomHelper.nextIntFromTo(0, 3) < 1){
			this.resident = 0; // Toggle whether fibroblasts are active based on necrosis
		}
	}
	
	// FIBROBLAST LEAVING, QUIESCENCE/HOMEOSTASIS
	public void fibrobQuiescence(){
		Context context = ContextUtils.getContext(this);
		int fibrobQuiTemp = (int)(fibroblastQuiescence);
		if (fibrobQuiTemp < 30){ // limits for fibroblast chance of quiescence/homeostasis
			fibrobQuiTemp = 30;
		}
		if (RandomHelper.nextIntFromTo(0, fibrobQuiTemp) < 1){ 
			// Either remove from the context OR make it a resident/stop moving and return to homeostatic/baseline conditions
			if (getResidentFibroblasts(context).size() < Math.ceil(Fiber.origFiberNumber*fibroScale*murphRep)){
				// If there is not a sufficient resident pool --> resident; otherwise chance to leave
				this.resident = 1;
			}
			else { // if there are enough residents apoptose/leave/migrate away
					context.remove(this);
			}
		}
	}
	
	// FIBROBLAST RECRUITMENT
	// Called from external source (Fiber) like SSC migration
	public static void fibrobRecruitment(Context<Object> context){
		List<Object> ecms = ECM.getECM(context); // Get a list of all the ecm agents to place cells on the ECM
			int fibRecTemp = (int)(70 - fibroblastRecruit/(Math.floor(Fiber.origFiberNumber))); // Fibroblast recruitment, weighted by recruitment signal
			if (fibRecTemp < 15){ 
				fibRecTemp = 15;
			}
			if (fibrobRecruitSaturation > getActiveFibroblasts(context).size() && RandomHelper.nextIntFromTo(0,fibRecTemp) < 1){ // if less than saturation point
				Fibroblast fibroblastNew = new Fibroblast(grid, space, 0, 0, 0, 1, 0, 0); // add a fibroblast
				context.add(fibroblastNew); // add to context
				// Get ECM location
				int index = RandomHelper.nextIntFromTo(0,  ecms.size() - 1); // draw random number with randomHelper
				Object ecmRandom = ecms.get(index); // get the ecm based on the random number chosen-- index within the list
				GridPoint ptECM = grid.getLocation(ecmRandom); // Get the ecm location
				grid.moveTo(fibroblastNew, ptECM.getX(), ptECM.getY()); // add the new fibroblast directly to the ECM
				((Fibroblast)fibroblastNew).setRecruitAge(1); // takes time to recruit				
			}
	}
	
	public void fibrobExpansion(){
		// Accounts for dependence of fibroblast expansion on satellite cells
		// Fibrob recruit occurs for approx first 48 hours, expansion lasts while recruit signal is high
		Context context = ContextUtils.getContext(this);
		List<Object> ecms = ECM.getECM(context); // Get a list of all the ecm agents to place cells on the ECM
		// Expansion is dependent on the number of satellite cells
		int normalizationNumTemp = Fiber.origFiberNumber;
		if (getFibroblasts(context).size() > Fiber.origFiberNumber){
			normalizationNumTemp = getFibroblasts(context).size(); // only use fibroblast number if greater than fiber number
		}
		int expChanceTemp = (int) (90 - fibroblastExpansion/normalizationNumTemp); // SSC expansion weighting factor- scaled by fibroblastExpansion signal // 83016 parameter SCALING
		if (expChanceTemp < 35){
			expChanceTemp = 35;
		}
		// fibroblast division chance-- assumes the fibroblasts divide to maintain counts-- more divisions = decreased chance of dividing- similar to SSC rules
		int divNumberEffect = 1;
			if (this.daughter >= 2){
				divNumberEffect = this.daughter + 1; // if the daughter has divided twice change the chance of dividing-- assume similar to other cells
			}
		if (this.expansionAge == 0 && fibroblastRecruit > Fiber.origFiberNumber*2 && (getFibroblasts(context).size() + getMyofibroblasts(context).size()) < fibroExpSaturation){
			// if under saturation limit
			// Expand with expansion signal (fibroblastExpansion); in the presence of recent damage (fibroblastRecruit); based on presence of SSC--lower dependence if normal amounts of fibroblasts
			// low chance of expansion with low signal:
			if (fibroblastExpansion < Fiber.origFiberNumber*10){ 
				int expChanceLow = (int)((200 - 10*fibroblastExpansion/Fiber.origFiberNumber));
				if (RandomHelper.nextIntFromTo(0, (int)(expChanceLow*Math.sqrt(divNumberEffect))) < 1){
					this.setExpansionAge(1); // this fibroblast starts proliferating- hyp: dividing like SSC
					this.setDaughter(this.getDaughter() + 1); // if the fibroblast divides, mark as a daughter
				}
			}
			// higher chance of expansion with high signal:
			else if (fibroblastExpansion > Fiber.origFiberNumber*10){
				if (RandomHelper.nextIntFromTo(0, (int)(expChanceTemp*Math.sqrt(divNumberEffect))) < 1){
					this.setExpansionAge(1); // this fibroblast starts proliferating- hyp: dividing like SSC
					this.setDaughter(this.getDaughter() + 1); // if the fibroblast divides, mark as a daughter
				}
			}
		}
		// If the fibroblast is dividing:
		if (this.expansionAge >= expTime){
			this.setExpansionAge(0);
			Fibroblast fibroblastNew = new Fibroblast(grid, space, 0, 0, 0, 0, 0, (this.getDaughter() + 1));  // add new fibroblast: mark as "daughter"
			context.add(fibroblastNew); // add to context
			// place next to current fibroblast
			// Get ECM location
			List<Object> neighbor = new ArrayList<Object>();
			VNQuery<Object> query = new VNQuery(grid, this, 1, 1); // get neighbors
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			// Choose a random neighbor and place the agent there
			for (Object neighborIter : iter){
				neighbor.add(neighborIter);
			}
			if (neighbor.size() > 0){
				int index = RandomHelper.nextIntFromTo(0,  neighbor.size() - 1); // draw random number with randomHelper
				Object randomNeigh = neighbor.get(index); // get the ecm based on the random number chosen-- index within the list
				GridPoint pt = grid.getLocation(randomNeigh); // Get the ecm location
				grid.moveTo(fibroblastNew, pt.getX(), pt.getY()); // add the new fibroblast directly to the ECM
			}
		}
	}
	
	// FIBROBLAST APOPTOSIS
	public void fibrobApoptosis(){
		Context context = ContextUtils.getContext(this);
		int origFiberNumber = Fiber.getOrigFiberNumber(); // Total number of muscle fibers, not elements with muscle in it
		// Chance of apoptosis vs. block apoptosis dependent on growth factors (signal above)
		if (fibroblastApop > fibroblastBlockApop*2 && fibroblastApop > 10*origFiberNumber && RandomHelper.nextIntFromTo(0,  (int)(55 - fibroblastApop/origFiberNumber)) < 2 && this.getApopAge() < fibrobApopCount){
			// Any cells- unless already apoptosing
			// if tnf is greater than tgf --> apoptosis of fibroblasts/nmmp/faps
			this.setApopAge(this.getApopAge() + 1); // will apop if it doesn't recruit first
		}
		else if (fibroblastBlockApop > fibroblastApop*2 && fibroblastBlockApop > 10*origFiberNumber && RandomHelper.nextIntFromTo(0,  (int)(120 - fibroblastBlockApop/origFiberNumber)) < 2 && this.getApopAge() < fibrobApopCount){ 
			// Any cells- unless already apoptosing
			// if tgf-beta is greatest, stop tnf-induced apoptosis of faps 
			setPhenotype(this.getPhenotype() + 1); // keep increasing phenotype number, at certain point --> myofibroblast
		}
		// If the two pressures are similar- either can happen:
		else if (fibroblastBlockApop/fibroblastApop < 2 && fibroblastBlockApop/fibroblastApop > .5 && RandomHelper.nextIntFromTo(0,  80) < 2 && this.getApopAge() < fibrobApopCount){
			if (RandomHelper.nextIntFromTo(0,1) < 1){
				this.setApopAge(this.getApopAge() + 1); // will apop if it doesn't recruit first
			}
			else {
				setPhenotype(this.getPhenotype() + 1); // keep increasing phenotype number, at certain point --> myofibroblast
			}
		}
		// APOPTOSIS
		if (this.getApopAge() > 12){	
			context.remove(this); // remove fibroblast
		}
	}
	
	public void fibrobNecrosisRemove(){
		Context context = ContextUtils.getContext(this);
		// Chance to remove random necrosis and replace with fibrosis after xx hours
		if (RandomHelper.nextIntFromTo(0, 9) < 1){
			List<Object> necrosis = Necrosis.getNecrosis(context);
			if (necrosis.size() != 0){
				int index = RandomHelper.nextIntFromTo(0,  necrosis.size() - 1);
				Object randomNecrosis = necrosis.get(index); // get the ecm based on the random number chosen-- index within the list
				GridPoint pt = grid.getLocation(randomNecrosis); // Get the ecm location
				context.remove(randomNecrosis); // remove necrosis
				// replace with lots of collagen
				ECM newECM = new ECM(grid, 10, 0, 0); // change to ECM
				context.add((ECM) newECM);
				grid.moveTo(newECM, pt.getX(), pt.getY()); 
			}
		}
	}
	
	// Myofibroblast secretions:
	public void myofibroblastSecretions(){
		// find location to add ecm from myofibroblast
		MooreQuery<Object> query = new MooreQuery(grid, this, 3, 3); // get neighbors in a wider area
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		// Go through the list of neighbors and find the neighbor with the lowest ECM
		Object ecmNeighbor = null;
		double tempLowColl = 5; // set a high collagen number to start- find values less than this
		for (Object neighbors : iter){
			if (neighbors instanceof ECM && ((ECM) neighbors).getCollagen() < tempLowColl){
				ecmNeighbor =neighbors;
				tempLowColl = ((ECM) neighbors).getCollagen();
			}
		}
		if (ecmNeighbor != null){
			((ECM) ecmNeighbor).setCollagen(((ECM) ecmNeighbor).getCollagen() + .2*(1/Fiber.tcf4Scale)); // add collagen to ecm element
		}
	}
	
	// Fibroblast move
	public void move(){ // Fibroblasts move towards necrosis/areas of low density collagen
		// At each step fibroblasts should sense surroundings and move toward necrosis/low collagen
		MooreQuery<Object> query = new MooreQuery(grid, this, 2, 2); // get neighbors in a wider area
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		List<Object> necrosisNeighbor = new ArrayList<Object>();
		List<Object> openNeighbor = new ArrayList<Object>();
		List<Object> lowCollNeighbor = new ArrayList<Object>();
		for (Object neighbor : iter){
			if (neighbor instanceof Necrosis || neighbor instanceof ECM){
				openNeighbor.add(neighbor); // if any neighbors are necrotic or ecm, add to the list
			}
			if (neighbor instanceof Necrosis){ // necrotic neighbors
				necrosisNeighbor.add(neighbor);
			}
			if (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() < 1){ // low collagen ecm neighbors
				lowCollNeighbor.add(neighbor);
			}
		}
		if (necrosisNeighbor.size() > 0 ){ // if there is necrosis move there
			int index = RandomHelper.nextIntFromTo(0, necrosisNeighbor.size()-1);
			Object randomNeighbor = necrosisNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		if (lowCollNeighbor.size() > 0 ){ // if there is low collagen
			int index = RandomHelper.nextIntFromTo(0, lowCollNeighbor.size()-1);
			Object randomNeighbor = lowCollNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
			// And secrete collagen
			if (((ECM) randomNeighbor).getCollagen() < .7*Fiber.mdxBaseCollagen){
				((ECM) randomNeighbor).setCollagen(((ECM) randomNeighbor).getCollagen() + .1*(1/Fiber.tcf4Scale));	
			}
			
		}
		else if(openNeighbor.size() > 0){ // Otherwise just pick a random direction of ecm go to it
			int index = RandomHelper.nextIntFromTo(0, openNeighbor.size()-1);
			Object randomNeighbor = openNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}		
	}
	
	public void moveTowards(GridPoint pt){
		GridPoint myPoint = grid.getLocation(this);
		GridPoint otherPoint = new GridPoint(pt.getX(), pt.getY());
		if ((pt.getY() - myPoint.getY()) > 1) {
			grid.moveTo(this, myPoint.getX(), myPoint.getY() + 1);			
		}
		else if ((pt.getY() - myPoint.getY()) < 1) {
			grid.moveTo(this, myPoint.getX(), myPoint.getY() - 1);
		}
		GridPoint myPointNew = grid.getLocation(this); // update location
		if ((pt.getX() - myPointNew.getX()) > 1) {
			grid.moveTo(this, myPointNew.getX() + 1, myPoint.getY());
		}
		else if ((pt.getX() - myPointNew.getX()) < 1) {
			grid.moveTo(this, myPointNew.getX() - 1, myPoint.getY());
		}
	}

	public static void initialize(Context<Object> context, Grid<Object> grid, int origFiberNumber, double fibrobMDX){
		// Start the model with a certain number of fibroblasts- written as an input value in the GUI
		for (int i = 0; i < Math.ceil(origFiberNumber*fibroScale*murphRep*fibrobMDX*Fiber.tcf4Scale); i++) {
			context.add(new Fibroblast(grid, space, 0, 0, 0, 1, 0, 0)); // Add the set number of fibroblasts to the context
		}
		List<Object> ecms = ECM.getECM(context); // Get a list of all the ecm agents to place cells on the ECM
		// Go through the lists of fibroblasts- at each fibroblast find a random ECM location and move the fibroblast there
		for (Object fibroblast : context){
			if (fibroblast instanceof Fibroblast){
				int index = RandomHelper.nextIntFromTo(0,  ecms.size() - 1); // draw random number with randomHelper
				Object ecmRandom = ecms.get(index); // get the ecm based on the random number chosen-- index within the list
				GridPoint ptECM = grid.getLocation(ecmRandom); // Get the ecm location
				grid.moveTo(fibroblast, ptECM.getX(), ptECM.getY());
			}
		}
	}
	
	public static List<Object> getFibroblasts(Context<Object> context){ // Get a list of all the fibroblasts
		List<Object> fibroblasts = new ArrayList<Object>(); // create a list of all the fibroblast agents
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).getPhenotype() < myofibroblastSwitch && ((Fibroblast) obj).getRecruitAge() == 0){
				// the recruitment time == transit time so the fibroblast is not technically on site yet
				fibroblasts.add(obj);
			}
		} 
		return fibroblasts;
	}
	
	public static List<Object> getMyofibroblasts(Context<Object> context){ // Get a list of all the myofibroblasts
		List<Object> myofibroblasts = new ArrayList<Object>(); // create a list of all the myofibroblast agents
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).phenotype >= myofibroblastSwitch){
				myofibroblasts.add(obj);
			}
		} 
		return myofibroblasts;
	}
	
	public static  List<Object> getActiveFibroblasts(Context<Object> context){ // Get a list of all the fibroblasts
		// Excludes "resident" fibroblasts and myofibroblasts
		List<Object> activefibroblasts = new ArrayList<Object>(); // create a list of all the fibroblast agents
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).resident == 0 && ((Fibroblast) obj).getRecruitAge() == 0 && ((Fibroblast) obj).phenotype < myofibroblastSwitch){ //
				activefibroblasts.add(obj);
			}
		} 
		return activefibroblasts;
	}
	
	
	public static  List<Object> getResidentFibroblasts(Context<Object> context){ // Get a list of all the fibroblasts
		// Get resident/quiescent fibrob population
		List<Object> resfibroblasts = new ArrayList<Object>(); // create a list of all the resident fibroblast agents
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).resident == 1 ){ //
				resfibroblasts.add(obj);
			}
		} 
		return resfibroblasts;
	}

	public double getFibroblastRecruit(){
		return fibroblastRecruit;
	}
	public double getFibroblastApop(){
		return fibroblastApop;
	}
	public double getFibroblastBlockApop(){
		return fibroblastBlockApop;
	}
	public double getFibroblastProlif(){
		return fibroblastProlif;
	}
	public double getFibroblastExpansion(){
		return fibroblastExpansion;
	}
	public double getFibroblastQuiescence(){
		return fibroblastQuiescence;
	}
	public static List<Object> getProlifFibroblasts(Context<Object> context){ 
		List<Object> prolifFibroblasts = new ArrayList<Object>(); 
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).recruitAge > 0){
				prolifFibroblasts.add(obj);
			}
		} 
		return prolifFibroblasts;
	}
	
	public double getProlifFibrobCount(){
		Context context = ContextUtils.getContext(this);
		return getProlifFibroblasts(context).size();
	}
	
	public static List<Object> getApopFibroblasts(Context<Object> context){ 
		List<Object> apopFibroblasts = new ArrayList<Object>(); 
		for (Object obj : context){ 
			if (obj instanceof Fibroblast && ((Fibroblast) obj).apopAge >= 2){
				apopFibroblasts.add(obj);
			}
		} 
		return apopFibroblasts;
	}
	
	public double getApopFibrobCount(){
		Context context = ContextUtils.getContext(this);
		return getApopFibroblasts(context).size();
	}
		
	public double eosinophilRecruit(){
		return eosinophilRecruit;
	}
	public void setResident(int resident){
		this.resident = resident;
	}
	public int getResident(){
		return resident;
	}	
	public void setPhenotype(int phenotype){
		this.phenotype = phenotype;
	}
	
	public int getPhenotype(){
		return phenotype;
	}
	public void setRecruitAge(int recruitAge){
		this.recruitAge = recruitAge;
	}
	
	public int getRecruitAge(){
		return recruitAge;
	}
	public void setApopAge(int apopAge){
		this.apopAge = apopAge;
	}
	
	public int getApopAge(){
		return apopAge;
	}
	public void setExpansionAge(int expansionAge){
		this.expansionAge = expansionAge;
	}
	
	public int getExpansionAge(){
		return expansionAge;
	}
	
	public double getFibroblastRecrSatTemp(){
		return fibrobRecruitSaturation;
	}
	
	public double getFibroExpSatTemp(){
		return fibroExpSaturation;
	}
	public void setDaughter(int daughter){
		this.daughter = daughter;
	}
	public int getDaughter(){
		return daughter;
	}
}

	
	
