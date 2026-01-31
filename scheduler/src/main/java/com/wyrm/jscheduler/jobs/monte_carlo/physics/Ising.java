package com.wyrm.jscheduler.jobs.monte_carlo.physics;

import com.wyrm.jscheduler.jobs.Job;
import org.apache.commons.math3.stat.descriptive.summary.Sum;
import java.util.ArrayList;
import java.util.Random;
import org.json.JSONObject;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadFactory;
import com.wyrm.jscheduler.utility.CustomFactory;
public class Ising implements Job, Runnable
{
    public static final  int INITIAL_SWEEPS = 1000;
    public static final  double BOLTZMAN_CONSTANT = 1;
    public static final  double MAX_SIZE = 256;

    private int[][] lattice;
    private int latticeSize;
    private final Random r = new Random();

    private JSONObject data;
    private int sweeps;
    private double temperature;
    private double absMagnetisation;
    private double coupling;
    private double externalMagnetisation;
    private ArrayList<Double> energySquared;
    private ArrayList<Double> energies;
    private ArrayList<Double> magnetisation;
    private ArrayList<Double> magnetisationSquared;
    public Ising()
    {}

    public Ising(JSONObject data)
    {
        this.data = data;
        int size = data.has("lattice_size") ? (Math.min(Integer.parseInt((String) data.get("lattice_size")), 256)) : 8;
        lattice = new int[size][size];
        if (((String)data.get("temperature_config")).equals("kelvin")) temperature = data.has("temperature") ? Double.parseDouble((String) data.get("temperature")):-273.13;
        else if (((String)data.get("temperature_config")).equalsIgnoreCase("celsius")) temperature = data.has("temperature") ? Double.parseDouble((String) data.get("temperature"))-273.13: -273.13;


        coupling = data.has("coupling") ? Double.parseDouble((String) data.get("coupling")):1.0;
        externalMagnetisation = data.has("magnetisation") ? Double.parseDouble((String) data.get("magnetisation")):0.0;
        sweeps = data.has("sweeps") ? Math.max(Math.abs(Integer.parseInt((String) data.get("sweeps"))),1000): 1000;
        energies = new ArrayList<>();
        magnetisation = new ArrayList<>();
        energySquared = new ArrayList<>();
        magnetisationSquared = new ArrayList<>();
        absMagnetisation = 0;
    }

    private void initialiseLattice()
    {
        for (int[] l: lattice)
        {
            for (int  i = 0 ;  i < l.length;i++)
            {
                l[i] = r.nextBoolean() ? 1:+1;

            }
        }
    }
    private void initialiseLattice(int value)
    {
        for (int[] l: lattice)
        {
            Arrays.fill(l, value);
        }
    }



    @Override
    public void run()
    {
        int count = 0;
        while (count < INITIAL_SWEEPS+sweeps)
        {
            metropolis(count);
            count++;
        }

    }

    /**
     * This method computes the sum of the neighbouring spin values of a specific point in the lattice
     * @param x The X co-ordinate, or simply, the 'row'.
     * @param y The Y co-ordinate, or simply, the 'column'.
     * @return Sum of neighbouring spin values
     */
    public double latticeNeighbours(int x , int y)
    {
        double up = (x-1 <0) ? lattice[latticeSize-1][y]: lattice[x-1][y];
        double down =  (x+1 >=latticeSize) ? lattice[0][y]: lattice[x+1][y];
        double left = (y-1 <0) ? lattice[x][lattice[x].length-1]: lattice[x][y-1];
        double right = (y+1 >= lattice[x].length) ? lattice[x][0] : lattice[x][y+1];
        return up+down+left+right;
    }
    /**
     * This method calculates the change in local energy for a specific spin value in the lattice
     * @param x The X co-ordinate, or simply, the 'row'.
     * @param y The Y co-ordinate, or simply, the 'column'.
     * @return Change in local energy
     */
    public double deltaEnergy(int x, int y)
    {
        return 2*coupling*-lattice[x][y]*(latticeNeighbours(x,y));
    }


    public JSONObject Monte()
    {
        JSONObject output;
        initialiseLattice();
        if (sweeps > 1000)
        {
            ExecutorService es = Executors.newFixedThreadPool((int) sweeps/1000);
        }
        int count = 0;
        while (count < INITIAL_SWEEPS+sweeps)
        {
            metropolis(count);
            count++;
        }


        double energyAverage = energies.stream().mapToDouble(Double::doubleValue).sum()/energies.size()/lattice.length*lattice.length;
        double energySquaredAverage = energySquared.stream().mapToDouble(Double::doubleValue).sum()/energies.size()*lattice.length*lattice.length;
        double magnetAverage = magnetisation.stream().mapToDouble(Double::doubleValue).sum()/energies.size()*(lattice.length*lattice.length);
        double magnetSquaredAverage= magnetisationSquared.stream().mapToDouble(Double::doubleValue).sum()/energies.size()*lattice.length*lattice.length;
        double heatCapacity = (energySquaredAverage - (energyAverage*energyAverage))/(BOLTZMAN_CONSTANT*temperature*temperature);
        double susceptibility = (magnetSquaredAverage - (magnetAverage*magnetAverage))/BOLTZMAN_CONSTANT*(temperature);

        output = new JSONObject();
        output.put("energy", energyAverage);
        output.put("energy_squared", energySquaredAverage);
        output.put("magnetisation", magnetAverage);
        output.put("magnetisation_squared", magnetSquaredAverage);
        output.put("absolute_mean_magnetisation", absMagnetisation/lattice.length/lattice.length);
        output.put("heat_capacity", heatCapacity);
        output.put("susceptibility", susceptibility);
        return output;
    }
    public void metropolis(int mcss)
    {
        for (int i = 0; i < lattice.length*lattice.length;i++)
        {
            int x = r.nextInt(0, lattice.length);
            int y = r.nextInt(lattice[0].length);
            double probability;
            double random;
            double energy = deltaEnergy(x, y);
            if (energy <= 0)
                lattice[x][y] = -lattice[x][y];
            else
            {
                probability = Math.exp(-energy /BOLTZMAN_CONSTANT * temperature);
                random = Math.random();
                synchronized(lattice)
                {
                    lattice[x][y] = (random < probability) ? -lattice[x][y]:lattice[x][y];
                }
            }
        }
        if (mcss > MIN_sweeps)
        {
            Sum s = new Sum();
            double magnetisation_ = 0;
            double hamiltonian = 0;
            for (int d = 0; d < lattice.length;d++)
            {
                magnetisation_ += Arrays.stream(lattice[d]).sum();
                for (int d_ = 0 ; d_ < lattice[0].length;d++)
                {
                    double down =  (x+1 >=latticeSize) ? lattice[0][y]: lattice[x+1][y];
                    double right = (y+1 >= lattice[x].length) ? lattice[x][0] : lattice[x][y+1];
                    hamiltonian+= lattice[d][d_] *(down+right);
                }
            }
            magnetisation_/=lattice.length*lattice.length;
            hamiltonian*=-coupling;
            magnetisation.add(magnetisation_);
            absMagnetisation += Math.abs(magnetisation_);
            energies.add(hamiltonian);
            energySquared.add(hamiltonian*hamiltonian);
            magnetisationSquared.add(magnetisation_*magnetisation_);
        }

    }
    public static void main(String[] args)
    {
//        double[][] aa = new double[3][3];
//        Arrays.stream(aa).forEach(row -> Arrays.fill(row, new Random().nextBoolean() ? 1:-1));
//        System.out.println(Arrays.deepToString(aa));

    }



}
