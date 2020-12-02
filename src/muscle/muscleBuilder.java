package muscle;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import com.jgoodies.common.base.Objects;

import au.com.bytecode.opencsv.CSVReader;
import repast.simphony.context.Context;
import repast.simphony.context.space.continuous.ContinuousSpaceFactory;
import repast.simphony.context.space.continuous.ContinuousSpaceFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.context.space.grid.GridFactory;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.RandomCartesianAdder;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.SimpleGridAdder;
import repast.simphony.space.grid.StrictBorders;

public class muscleBuilder implements ContextBuilder<Object> {

	// Context builder builds a context. Context is a named set of agents. Building one = naming it and adding agents.
	// Projections take the agents in a Context and impose some sort of structure on them, such as continuousSpace and Grid.
	// This context will hold the fibroblasts, muscle cells, ecm, and satellite cells
	
	public static final int gridx = 150; // x dimension of abm grid
	public static final int gridy = 190; // y dimension of abm grid
	
	@Override
	public Context<Object> build(Context<Object> context) {
		context.setId("muscle");
		// Create grid and continuous space for muscle context
		
		// Note: Grid size should be larger than the imported geometry
		GridFactory gridFactory = GridFactoryFinder.createGridFactory(null);
		Grid<Object> grid = gridFactory.createGrid("grid", context, new GridBuilderParameters<Object>(new StrictBorders(), new SimpleGridAdder<Object>(), true, gridx, gridy));		
		
		ContinuousSpaceFactory spaceFactory = ContinuousSpaceFactoryFinder.createContinuousSpaceFactory(null);
		ContinuousSpace<Object> space = spaceFactory.createContinuousSpace("space", context, new RandomCartesianAdder<Object>(), 
				new repast.simphony.space.continuous.StrictBorders(), gridx, gridy); 
		
		// Create a network
		NetworkBuilder<Object> netBuilder = new NetworkBuilder<Object>("fiber network", context, true);
		netBuilder.buildNetwork(); // use network builder to create the network, call buildnetwork to actually build the network
		
		// Read in a CSV file for the fiber and ecm agents: CSV File: xCoor, yCoor, Material Type (1 = border, 2 = ecm, 3 = fiber)
				try {
			readLineByLine(context, grid);
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Add a single "InflamCell" to the context for solving
		for (int i = 0; i < 1; i++) {
			context.add(new InflamCell(grid, space)); 
		}
		// Add a single "GrowthFactor" cell to the context for solving and printing data
		for (int i = 0; i < 1; i++) {
			context.add(new GrowthFactors(grid, space)); 
		}
		
		// Set duration of simulation:
		RunEnvironment.getInstance().endAt(672);
		//RunEnvironment.getInstance().endAt(2160);
		
		return context; // return the muscle context
		
		
	}
	
	public static void readLineByLine(Context<Object> context, Grid<Object> grid) throws IOException {
		// Read in CSV file and convert it to the initial geometry read by the ABM
		// if disease state == 0 --> healthy
		String csvFilename = "./data/healthy_histo.csv";
		// If disease state- rename
		if (Fiber.diseaseState == 1){
			csvFilename = "./data/youngMDX_histo.csv";
		}
		else if (Fiber.diseaseState == 2){
			csvFilename = "./data/adultMDX_histo.csv";
		}
		else if (Fiber.diseaseState == 3){
			csvFilename = "./data/oldMDX_histo.csv";
		}
		
		CSVReader csvReader = new CSVReader(new FileReader(csvFilename));
		String[] row = null;

		while ((row = csvReader.readNext()) != null) {
			if ((row[2].equals("3"))) { // type 3 = muscle; if muscle add a fiber
				Fiber fiber = new Fiber(grid, 0, 0, 0, 0, 0, 0) ;
				context.add(fiber) ;
				grid.moveTo(fiber, Integer.parseInt(row[0]), Integer.parseInt(row[1]));
			}
			//if (Objects.equals(row[2], "1")) { // CTO paper
			if ((row[2].equals("2"))) { // type 2 = ecm; if muscle add a fiber
				ECM ecm = new ECM(grid, 1.0, 0, 0) ;
				context.add(ecm) ;
				grid.moveTo(ecm, Integer.parseInt(row[0]), Integer.parseInt(row[1]));
			}
		}
		csvReader.close();
	} 
}



			




