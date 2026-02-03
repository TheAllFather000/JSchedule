package com.wyrm.jscheduler.jobs.monte_carlo.physics;

import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.json.JSONObject;

import com.wyrm.jscheduler.jobs.Job;
import com.wyrm.jscheduler.utility.CustomFactory;
public class Ising implements Job, Callable<Double[]>
{
    public static final  int INITIAL_SWEEPS = 5000;
    public static final  int INITIAL_SWEEPS_LARGER = 10000;
    public static final  double BOLTZMAN_CONSTANT = 1;
    public static final  int MAX_SIZE = 256;

    private int[][] lattice;
    private int latticeSize;
    private final Random r = new Random();

    private JSONObject data;
    private int sweeps;
    private int measure;
    private double temperature;
    private double absMagnetisation;
    private double coupling;
    private double binderCumulant;
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
        latticeSize = data.has("lattice_size") ? (Math.min(Integer.parseInt((String) data.get("lattice_size")), MAX_SIZE)) : 8;
        lattice = new int[latticeSize][latticeSize];
        temperature = data.has("temperature") ? Double.parseDouble((String) data.get("temperature")):1.5;
        temperature = Math.abs(temperature);

        coupling = data.has("coupling") ? Double.parseDouble((String) data.get("coupling")):1.0;
        coupling = Math.abs(coupling);
//        externalMagnetisation = data.has("magnetisation") ? Double.parseDouble((String) data.get("magnetisation")):0.0;
        sweeps = data.has("sweeps") ? Math.max(Math.abs(Integer.parseInt((String) data.get("sweeps"))),1000): 1000;
        measure = data.has("measure") ? Math.abs(Integer.parseInt((String)data.get("measure"))):10;
        if (measure ==0)
            measure = 10;
        energies = new ArrayList<>();
        magnetisation = new ArrayList<>();
        energySquared = new ArrayList<>();
        magnetisationSquared = new ArrayList<>();
        absMagnetisation = 0;
        binderCumulant = 0;
    }
    public Ising(JSONObject data, int[][] l)
    {
        this.data = data;
        latticeSize = data.has("lattice_size") ? (Math.min(Integer.parseInt((String) data.get("lattice_size")), MAX_SIZE)) : 8;
        lattice = l;
        temperature = data.has("temperature") ? Double.parseDouble((String) data.get("temperature")):1.5;
        temperature = Math.abs(temperature);

        coupling = data.has("coupling") ? Double.parseDouble((String) data.get("coupling")):1.0;
        coupling = Math.abs(coupling);
//        externalMagnetisation = data.has("magnetisation") ? Double.parseDouble((String) data.get("magnetisation")):0.0;
        sweeps = data.has("sweeps") ? Math.max(Math.abs(Integer.parseInt((String) data.get("sweeps"))),1000): 1000;
        measure = data.has("measure") ? Math.abs(Integer.parseInt((String)data.get("measure"))):10;
        if (measure ==0)
            measure = 10;
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
                l[i] = r.nextBoolean() ? 1:-1;

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
    public Double[] call()
    {
        initialiseLattice();
        int initialSweeps = (latticeSize >= 32) ? INITIAL_SWEEPS_LARGER : INITIAL_SWEEPS;
        int totalSweepsPerThread = initialSweeps + (sweeps / Runtime.getRuntime().availableProcessors());
        System.out.println(initialSweeps + (sweeps / Runtime.getRuntime().availableProcessors()));
        int count = 0;
        while (count < totalSweepsPerThread)
        {
            metropolis(count, initialSweeps);
            count++;
        }
        Double[] d = new Double[5];

        int N = lattice[0].length * lattice.length;
        double energySum = energies.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double energySquaredSum = energySquared.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double magnetSum = magnetisation.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double magnetSquaredSum = magnetisationSquared.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        d[0] = energySum;
        d[1] = energySquaredSum;
        d[2] = magnetSum;
        d[3] = magnetSquaredSum;
        d[4] = absMagnetisation / magnetisation.size();
        return d;
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
        return 2*coupling*lattice[x][y]*(latticeNeighbours(x,y));
    }


    public JSONObject Monte() throws ExecutionException, InterruptedException {
        JSONObject output = null;
        if (latticeSize >= 100 && sweeps >1000)
        {
            ArrayList<Callable<Double[]>> todo = new ArrayList<>();
            for (int i = 0; i < Runtime.getRuntime().availableProcessors();i++)
                todo.add(new Ising(data));
            System.out.println("balls");
            ExecutorService es = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors(), new CustomFactory());
            ArrayList<Future<Double[]>> results = (ArrayList<Future<Double[]>>) es.invokeAll(todo);
            try {
                es.shutdown();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            double energy = 0;
            double energySquared = 0;
            double magnet = 0;
            double magnetSquared = 0;
            double absMagnetisation = 0;
            for (Future<Double[]> result : results)
            {
                try {
                    Double[] f = result.get();
                    energy += f[0];
                    energySquared += f[1];
                    magnet += f[2];
                    magnetSquared += f[3];
                    absMagnetisation += f[4];
                } catch (Exception e) {
                    e.printStackTrace();
                }


            }
            energy /= results.size();
            magnet /= results.size();
            magnetSquared /= results.size();
            absMagnetisation /= results.size();
            energySquared /=results.size();


            int N = lattice.length*lattice[0].length;
            output = new JSONObject();
            output = new JSONObject();
            output.put("energy", energy);
            output.put("energy_per_spin", energy/(lattice.length*lattice[0].length));
            output.put("energy_squared", energySquared);
            output.put("magnetisation", magnet);
            output.put("magnetisation_per_spin", Math.abs(magnet)/(lattice.length*lattice[0].length));
            output.put("magnetisation_squared", magnetSquared);
            output.put("absolute_mean_magnetisation", absMagnetisation);
            output.put("heat_capacity", (energySquared- (energy*energy))/(BOLTZMAN_CONSTANT*temperature*temperature*N));
            output.put("susceptibility", (magnetSquared-(magnet*magnet))/(BOLTZMAN_CONSTANT*temperature*N));
            System.out.println("GGG");
            return output;

        }
        else {
            initialiseLattice(1);
            int count = 0;
            int initialSweeps = (latticeSize >= 32) ? INITIAL_SWEEPS_LARGER : INITIAL_SWEEPS;
            while (count < initialSweeps + sweeps) {
                metropolis(count, initialSweeps);
                count++;
            }

            int N = lattice.length*lattice[0].length;
            double energyAverage = energies.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double energySquaredAverage = energySquared.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double magnetAverage = magnetisation.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double magnetSquaredAverage = magnetisationSquared.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double heatCapacity = (energySquaredAverage - (energyAverage * energyAverage)) / (BOLTZMAN_CONSTANT * temperature * temperature*N);
            double susceptibility = (magnetSquaredAverage - (magnetAverage * magnetAverage)) / (BOLTZMAN_CONSTANT * temperature*N);

            output = new JSONObject();
            binderCumulant = 1- (binderCumulant/(3*magnetSquaredAverage*magnetSquaredAverage));
            output.put("energy", energyAverage);
            output.put("energy_per_spin", energyAverage/(N));
            output.put("energy_squared", energySquaredAverage);
            output.put("magnetisation", magnetAverage);
            output.put("magnetisation_per_spin", magnetAverage/(N));
            output.put("binder_cumulant", binderCumulant);
            output.put("magnetisation_squared", magnetSquaredAverage);
            System.out.println(energyAverage);
            System.out.println(energySquaredAverage);
            output.put("absolute_mean_magnetisation", absMagnetisation / magnetisation.size());
            output.put("heat_capacity", heatCapacity);
            output.put("susceptibility", susceptibility);
            return output;
        }
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
                probability = Math.exp(-energy /(BOLTZMAN_CONSTANT * temperature));
                random = r.nextDouble();
                lattice[x][y] = (random < probability) ? -lattice[x][y]:lattice[x][y];
            }
        }
        if (mcss > INITIAL_SWEEPS && mcss%measure ==0)
        {
            double magnetisation_ = 0;
            double hamiltonian = 0;
            for (int d = 0; d < lattice.length;d++)
            {
                magnetisation_ += Arrays.stream(lattice[d]).sum();
                for (int d_ = 0 ; d_ < lattice[0].length;d_++)
                {
                    double down =  (d+1 >=latticeSize) ? lattice[0][d_]: lattice[d+1][d_];
                    hamiltonian += -coupling * lattice[d][d_] * down;
                    double right = (d_+1 >= lattice[d].length) ? lattice[d][0] : lattice[d][d_+1];
                    hamiltonian += -coupling * lattice[d][d_] * right;

                }
            }
            binderCumulant +=Math.pow(magnetisation_, 4);
            magnetisation.add(magnetisation_);
            absMagnetisation += Math.abs(magnetisation_);
            energies.add(hamiltonian);
            energySquared.add(hamiltonian*hamiltonian);
            magnetisationSquared.add(magnetisation_*magnetisation_);
        }

    }
    public void metropolis(int mcss, int sweep)
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
                probability = Math.exp(-energy /(BOLTZMAN_CONSTANT * temperature));
                random = r.nextDouble();
                lattice[x][y] = (random < probability) ? -lattice[x][y]:lattice[x][y];
            }
        }
        if (mcss > sweep && mcss%measure ==0)
        {
            double magnetisation_ = 0;
            double hamiltonian = 0;
            for (int d = 0; d < lattice.length;d++)
            {
                magnetisation_ += Arrays.stream(lattice[d]).sum();
                for (int d_ = 0 ; d_ < lattice[0].length;d_++)
                {
                    double down =  (d+1 >=latticeSize) ? lattice[0][d_]: lattice[d+1][d_];
                    hamiltonian += -coupling * lattice[d][d_] * down;
                    double right = (d_+1 >= lattice[d].length) ? lattice[d][0] : lattice[d][d_+1];
                    hamiltonian += -coupling * lattice[d][d_] * right;

                }
            }
            binderCumulant +=Math.pow(magnetisation_, 4);
            magnetisation.add(magnetisation_);
            absMagnetisation += Math.abs(magnetisation_);
            energies.add(hamiltonian);
            energySquared.add(hamiltonian*hamiltonian);
            magnetisationSquared.add(magnetisation_*magnetisation_);
        }

    }
    //    public void metropolis(int mcss, int[][] l)
//    {
//        for (int i = 0; i < l.length*l[0].length;i++)
//        {
//            int x = r.nextInt(0, l.length);
//            int y = r.nextInt(l[0].length);
//            double probability;
//            double random;
//            double energy = deltaEnergy(x, y);
//            if (energy <= 0)
//                l[x][y] = -l[x][y];
//            else
//            {
//                probability = Math.exp(-energy /(BOLTZMAN_CONSTANT * temperature));
//                random = Math.random();
//                l[x][y] = (random < probability) ? -l[x][y]:l[x][y];
//            }
//        }
//        if (mcss > INITIAL_SWEEPS && mcss%measure ==0)
//        {
//            double magnetisation_ = 0;
//            double hamiltonian = 0;
//            for (int d = 0; d < l.length;d++)
//            {
//                magnetisation_ += Arrays.stream(l[d]).sum();
//                for (int d_ = 0 ; d_ < l[0].length;d_++)
//                {
//                    double down =  (d+1 >=l.length) ? l[0][d_]: l[d+1][d_];
//                    hamiltonian += -coupling * lattice[d][d_] * down;
//                    double right = (d_+1 >= l[d].length) ? l[d][0] : l[d][d_+1];
//                    hamiltonian += -coupling * lattice[d][d_] * right;
//                }
//            }
//            hamiltonian*=-coupling;
//            System.out.println("Hami "+hamiltonian);
//            magnetisation.add(magnetisation_);
//            absMagnetisation += Math.abs(magnetisation_);
//            energies.add(hamiltonian);
//            energySquared.add(hamiltonian*hamiltonian);
//            magnetisationSquared.add(magnetisation_*magnetisation_);
//        }
//
//    }
    public static void main(String[] args)
    {
//        double[][] aa = new double[3][3];
//        Arrays.stream(aa).forEach(row -> Arrays.fill(row, new Random().nextBoolean() ? 1:-1));
//        System.out.println(Arrays.deepToString(aa));
        JSONObject o = new JSONObject();
        o.put("temperature", "5.0");
        o.put("sweeps", "50000");
        o.put("lattice_size", "100");
        o.put("coupling", "1");
        o.put("measure", "100");
        System.out.println(Runtime.getRuntime().availableProcessors());
        Ising i = new Ising(o);
        Instant s = Instant.now();

        try
        {
            o = i.Monte();
        }
        catch (Exception e)
        {
            e.printStackTrace();;
        }
        System.out.println(o);
        Instant f =Instant.now();
        System.out.println(Duration.between(s, f).toSeconds());

    }



}
