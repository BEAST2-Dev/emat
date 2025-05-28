package test.emat.stochasticmapping;

import java.util.List;

import org.junit.Test;

import beast.base.util.Randomizer;
import emat.stochasticmapping.TimeStateInterval;
import emat.stochasticmapping.UniformisationStochasticMapping;

public class TestUniformizationStochasticMapping {

	@Test
    public void test() {
        // Example from Wikipedia for CTMC
        // States: 0 (E), 1 (A)
        // R = [[-alpha, alpha], [beta, -beta]]
        double alpha = 0.07;
        double beta = 0.03;
        double[][] R = {
                {-alpha, alpha},
                {beta, -beta}
        };
        
        // Calculate M_MAX_JUMPS dynamically (example)
        double lambdaMaxCalc = Math.max(-R[0][0], -R[1][1]);
        double totalTimeExample = 10.0;
        int mMaxJumpsDynamic = (int) Math.ceil(lambdaMaxCalc * totalTimeExample + 10 * Math.sqrt(lambdaMaxCalc * totalTimeExample));
        if (mMaxJumpsDynamic < 20) mMaxJumpsDynamic = 20; // a minimum
        if (mMaxJumpsDynamic > 20) mMaxJumpsDynamic = 20; // a practical cap to avoid huge factorial/power issues

        System.out.println("Dynamic M_MAX_JUMPS: " + mMaxJumpsDynamic);

        Randomizer.setSeed(1271);
        UniformisationStochasticMapping usm = new UniformisationStochasticMapping(mMaxJumpsDynamic);

        int startState = 0;
        int endState = 1; // Try 0 -> 0 and 0 -> 1

        System.out.println("\n--- Mapping from " + startState + " to " + endState + " in time " + totalTimeExample + " ---");
        List<TimeStateInterval> path1 = usm.generatePath(R, startState, endState, totalTimeExample);
        if (path1 != null) {
            path1.forEach(System.out::println);
            double totalDuration = path1.stream().mapToDouble(p -> p.endTime() - p.startTime()).sum();
            System.out.println("Total path duration: " + totalDuration);
        } else {
            System.out.println("Path could not be generated.");
        }

        System.out.println("\n--- Mapping from " + startState + " to " + startState + " in time " + totalTimeExample + " ---");
        List<TimeStateInterval> path2 = usm.generatePath(R, startState, startState, totalTimeExample);
         if (path2 != null) {
            path2.forEach(System.out::println);
            double totalDuration = path2.stream().mapToDouble(p -> p.endTime() - p.startTime()).sum();
            System.out.println("Total path duration: " + totalDuration);
        } else {
            System.out.println("Path could not be generated.");
        }


        // More complex example (3 states)
        double[][] R3 = {
            {-0.3,  0.2,  0.1},
            { 0.1, -0.5,  0.4},
            { 0.05, 0.15, -0.2}
        };
        double totalTime3 = 10.0;
        startState = 0;
        endState = 2;
        
        lambdaMaxCalc = 0;
        for(int i=0; i<R3.length; ++i) lambdaMaxCalc = Math.max(lambdaMaxCalc, -R3[i][i]);
        mMaxJumpsDynamic = (int) Math.ceil(lambdaMaxCalc * totalTime3 + 10 * Math.sqrt(lambdaMaxCalc * totalTime3));
        if (mMaxJumpsDynamic < 20) mMaxJumpsDynamic = 20; 
        if (mMaxJumpsDynamic > 20) mMaxJumpsDynamic = 20;

        UniformisationStochasticMapping usm3 = new UniformisationStochasticMapping(mMaxJumpsDynamic);

        System.out.println("\n--- 3-State: Mapping from " + startState + " to " + endState + " in time " + totalTime3 + " ---");
        List<TimeStateInterval> path3 = usm3.generatePath(R3, startState, endState, totalTime3);
        if (path3 != null) {
            path3.forEach(System.out::println);
             double totalDuration = path3.stream().mapToDouble(p -> p.endTime() - p.startTime()).sum();
            System.out.println("Total path duration: " + totalDuration);
        } else {
            System.out.println("Path could not be generated.");
        }
    }
}
