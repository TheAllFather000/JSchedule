package com.wyrm.jscheduler.jobs.monte_carlo.physics;
import java.time.Duration;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import lombok.Data;
import org.json.JSONObject;
import com.wyrm.jscheduler.jobs.Job;
import com.wyrm.jscheduler.utility.CustomFactory;
@Data

public class Ising implements Job, Callable<Double[]>
{
    public static final  int INITIAL_SWEEPS = 5000;
    public static final  int INITIAL_SWEEPS_LARGER = 10000;
    public static final  double BOLTZMAN_CONSTANT = 1;
    public static final  int MAX_SIZE = 256;
    public static int THERM = 10;
    public static int THERM_LARGER = 100;
    public static int MEASURES = 10000;
    public static int MEASURES_LARGER = 100000;
    public static final int[] LATTICE_SIZES = new int[]{16, 32, 64, 128, 256};

    private int[][] lattice;
    private int latticeSize;
    private final Random r = new Random();
    private JSONObject data;
    private int sweeps;
    private int row;
    private int col;
    private int measure;
    private double temperature;
    private double absMagnetisation;
    private double coupling;
    private double binderCumulant;
    private Stack<int[]> stack;
    private double externalMagnetisation;
    private boolean[][] cluster;
    private ArrayList<Double> energies;
    private ArrayList<Double> energySquared;
    private ArrayList<Double> magnetisation;
    private ArrayList<Double> magnetisationSquared;
    private ArrayList<Double> magnetisation4;
    private ArrayList<Double> magnetisation2;
    public Ising()
    {}
    public Ising(JSONObject data)
    {
        this.data = data;
        row = data.has("row") ? Math.abs(Integer.parseInt((String) data.get("row"))) : 8;
        col = data.has("col") ? Math.abs(Integer.parseInt((String) data.get("col"))) : 8;
        lattice = new int[row][col];
        stack = new Stack<>();
        temperature = data.has("temperature") ? Double.parseDouble((String) data.get("temperature")):1.5;
        temperature = Math.abs(temperature);
        coupling = data.has("coupling") ? Double.parseDouble((String) data.get("coupling")):1.0;
        coupling = Math.abs(coupling);
        cluster = new boolean[row][col];
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
        magnetisation4 = new ArrayList<>();
        magnetisation2 = new ArrayList<>();
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
    private void initialiseLattice(int row, int col)
    {
        this.row = row;
        this.col = col;
        lattice = new int[row][col];
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
        int initialSweeps = (row*col >= 1000) ? INITIAL_SWEEPS_LARGER : INITIAL_SWEEPS;
        int totalSweepsPerThread = initialSweeps + (sweeps / Runtime.getRuntime().availableProcessors());
        System.out.println(initialSweeps + (sweeps / Runtime.getRuntime().availableProcessors()));
        int count = 0;
        while (count < totalSweepsPerThread)
        {
            metropolis(count, initialSweeps);
            count++;
        }
        Double[] d = new Double[7];

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
        d[5] = magnetisation2.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        d[6] = magnetisation4.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        return d;
    }




    /**
     * This method computes the sum of the neighbouring spin values of a specific point in the lattice
     * @param x The X co-ordinate, or simply, the 'row'.
     * @param y The Y co-ordinate, or simply, the 'column'.
     * @return Sum of neighbouring spin values
     */
    public int latticeNeighbours(int x , int y)
    {
        int up = (x-1 <0) ? lattice[row-1][y]: lattice[x-1][y];
        int down =  (x+1 >=lattice.length) ? lattice[0][y]: lattice[x+1][y];
        int left = (y-1 <0) ? lattice[x][lattice[x].length-1]: lattice[x][y-1];
        int right = (y+1 >= lattice[x].length) ? lattice[x][0] : lattice[x][y+1];
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


    public JSONObject Metro() throws ExecutionException, InterruptedException {
        JSONObject output = null;
        if (row*col >= 10000 && sweeps >1000)
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
            double m4 = 0;
            double m2 = 0;
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
                    m2 +=f[5];
                    m4 += f[6];
                } catch (Exception e) {
                    e.printStackTrace();
                }


            }
            energy /= results.size();
            magnet /= results.size();
            magnetSquared /= results.size();
            absMagnetisation /= results.size();
            energySquared /=results.size();
            m2 /=results.size();
            m4/=results.size();

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
            output.put("binder_cumulant", 1- ((m4)/(3*m2*m2)));
            System.out.println("GGG");
            return output;

        }
        else {
            initialiseLattice();
            int count = 0;
            int initialSweeps = (row*col >= 1000) ? INITIAL_SWEEPS_LARGER : INITIAL_SWEEPS;
            while (count < initialSweeps + sweeps) {
                metropolis(count, initialSweeps);
                count++;
            }

            int N = lattice.length*lattice[0].length;
            double energyAverage = energies.parallelStream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double energySquaredAverage = energySquared.parallelStream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double magnetAverage = magnetisation.parallelStream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double magnetSquaredAverage = magnetisationSquared.parallelStream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double heatCapacity = (energySquaredAverage - (energyAverage * energyAverage)) / (BOLTZMAN_CONSTANT * temperature * temperature*N);
            double susceptibility = (magnetSquaredAverage - (magnetAverage * magnetAverage)) / (BOLTZMAN_CONSTANT * temperature*N);
            double m2 = magnetisation2.parallelStream().mapToDouble(Double::doubleValue).average().orElse(0);
            double m4 = magnetisation4.parallelStream().mapToDouble(Double::doubleValue).average().orElse(0);

            output = new JSONObject();

            binderCumulant = 1- ((m4)/(3*m2*m2));
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
        for (int i = 0; i < lattice.length*lattice[0].length;i++)
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
                    double down =  (d+1 >=row) ? lattice[0][d_]: lattice[d+1][d_];
                    hamiltonian += -coupling * lattice[d][d_] * down;
                    double right = (d_+1 >= lattice[d].length) ? lattice[d][0] : lattice[d][d_+1];
                    hamiltonian += -coupling * lattice[d][d_] * right;

                }
            }
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
                    double down =  (d+1 >=row) ? lattice[0][d_]: lattice[d+1][d_];
                    hamiltonian += -coupling * lattice[d][d_] * down;
                    double right = (d_+1 >= lattice[d].length) ? lattice[d][0] : lattice[d][d_+1];
                    hamiltonian += -coupling * lattice[d][d_] * right;

                }
            }
            magnetisation.add(magnetisation_);
            absMagnetisation += Math.abs(magnetisation_);
            energies.add(hamiltonian);
            energySquared.add(hamiltonian*hamiltonian);
            magnetisationSquared.add(magnetisation_*magnetisation_);
            int N = (row*col);
            double m = magnetisation_/N;
            magnetisation4.add(m*m*m*m);
            magnetisation2.add(m*m);
        }

    }

    public JSONObject wolff()
    {
        int count = 0;
        Thread t;
        while (count <= sweeps)
        {
            cluster();
            if (count >=Math.max(10, Math.log(col*row)))
                measure();
            count++;
        }
        int N = lattice.length*lattice[0].length;
        JSONObject output = new JSONObject();
        double energy = energies.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double energySquare = energySquared.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double magnet = magnetisation.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double magnet2 = magnetisationSquared.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double m4 = magnetisation4.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double m2 = magnetisation2.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double heatCapacity = (energySquare - (energy * energy)) / (BOLTZMAN_CONSTANT * temperature * temperature*N);
        double susceptibility = (magnet2 - (magnet * magnet)) / (BOLTZMAN_CONSTANT * temperature*N);

        binderCumulant = 1-(m4/(3*m2*m2));
//        System.out.println(energy);
//        System.out.println(energySquare);
//        System.out.println(magnet);
//        System.out.println(magnet2);
//        System.out.println(heatCapacity);
//        System.out.println(susceptibility);
//        System.out.println(binderCumulant);
        output.put("energy", energy);
        output.put("energy_per_spin", energy/(N));
        output.put("energy_squared", energySquare);
        output.put("magnetisation", magnet);
        output.put("magnetisation_per_spin", magnet/(N));
        output.put("binder_cumulant", binderCumulant);
        output.put("magnetisation_squared", magnet2);
        output.put("absolute_mean_magnetisation", absMagnetisation / magnetisation.size());
        output.put("heat_capacity", heatCapacity);
        output.put("susceptibility", susceptibility);
        return output;
               
    }
    public void wolff(boolean t) {
        if (row > 64) {
            for (int i = 0; i < THERM; i++)
                cluster();
            for (int i = 0; i < MEASURES; i++) {
                cluster();
                if (i % row == 0)
                    measure(true);
            }

        }
        else
        {
            for (int i = 0; i < THERM_LARGER; i++)
                cluster();
            for (int i = 0; i < MEASURES_LARGER; i++) {
                cluster();
                if (i % row == 0)
                    measure(true);
            }
        }
        double m4 = magnetisation4.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double m2 = magnetisation2.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        binderCumulant = 1 - (m4 / (3 * m2 * m2));

    }
    public void measure()
    {
        double magnetisation_ = 0;
        double hamiltonian = 0;
        for (int d = 0; d < lattice.length;d++)
        {
            magnetisation_ += Arrays.stream(lattice[d]).sum();
            for (int d_ = 0 ; d_ < lattice[0].length;d_++)
            {
                double down =  (d+1 >=lattice.length) ? lattice[0][d_]: lattice[d+1][d_];
                hamiltonian += -coupling * lattice[d][d_] * down;
                double right = (d_+1 >= lattice[d].length) ? lattice[d][0] : lattice[d][d_+1];
                hamiltonian += -coupling * lattice[d][d_] * right;

            }
        }
        absMagnetisation += Math.abs(magnetisation_);
        energies.add(hamiltonian);
        magnetisation.add(magnetisation_);
        energySquared.add(hamiltonian*hamiltonian);
        magnetisationSquared.add(magnetisation_*magnetisation_);
        int N = (row*col);
        double m = magnetisation_/N;
        magnetisation4.add(m*m*m*m);
        magnetisation2.add(m*m);
    }

    public void measure(boolean t)
    {
        double magnetisation_ = 0;
        double hamiltonian = 0;
        for (int[] d : lattice)
        {
            magnetisation_ += Arrays.stream(d).sum();
        }
        int N = (row*col);
        double m = magnetisation_/N;
        magnetisation4.add(m*m*m*m);
        magnetisation2.add(m*m);
    }
    public int[][] getSiteNeighbours(int x, int y)
    {
        int[][] i = new int[4][2];
//        int up = (x-1 <0) ? lattice[row-1][y]: lattice[x-1][y];
//        int down =  (x+1 >=row) ? lattice[0][y]: lattice[x+1][y];
//        int left = (y-1 <0) ? lattice[x][lattice[x].length-1]: lattice[x][y-1];
//        int right = (y+1 >= lattice[x].length) ? lattice[x][0] : lattice[x][y+1];
        i[0][0] = (x-1 < 0) ? row-1:x-1; i[0][1] = y;
        i[1][0] = (x+1 >=row) ? 0:x+1; i[1][1] = y;
        i[2][0] = x; i[2][1] = (y-1 < 0) ? col-1:y-1;
        i[3][0] = x; i[3][1] = (y+1 >= col) ? 0: y+1;

        return i;
    }
    public void cluster()
    {
        int x = r.nextInt(row);
        int y = r.nextInt(col);
        stack.push(new int[]{x, y});
        cluster[x][y] = true;
        int spinVal = lattice[x][y];
        double probability = 1.0 - Math.exp(-2*(coupling/temperature));
        while (!stack.isEmpty())
        {
            int[] site = stack.pop();
            int[][] neighbours = getSiteNeighbours(site[0], site[1]);
            for (int[] i : neighbours)
            {
                if (lattice[i[0]][i[1]] == spinVal && !cluster[i[0]][i[1]])
                {
                    double rand = r.nextDouble();
                    if (rand < probability)
                    {
                        cluster[i[0]][i[1]] = true;
                        stack.push(i);
                    }
                }
            }

        }
        for (int i = 0; i < cluster.length;i++)
        {
            for (int j = 0; j < cluster[0].length;j++)
            {
                if (cluster[i][j])
                {
                    lattice[i][j]*=-1;
                }
            }
        }

        cluster = new boolean[row][col];
    }

    public HashMap<Integer, String> binderCrossing()  throws InterruptedException, RuntimeException
    {
        JSONObject output = new JSONObject();
        HashMap<Integer, String> Datapoints = new HashMap<>();
        int index = 0;
        double o = 2.0;
        if (sweeps < 1000)
        {
            for (int L : LATTICE_SIZES)
            {
                Ising ing = new Ising();
                initialiseLattice(L, L);
                double count = 2.0;
                while (count < 2.6)
                {
                    row = L;
                    col = L;
                    temperature = count;
                    ing.magnetisation4 = new ArrayList<>();
                    ing.magnetisation2 = new ArrayList<>();
                    wolff(true);
                    if (Datapoints.containsKey(L)) Datapoints.put(L,Datapoints.get(L) + ", "+Arrays.toString(new Double[]{count, binderCumulant}));
                    else Datapoints.put(L, Arrays.toString(new Double[]{count, binderCumulant}));
                    count+=0.01;
                }
            }
        }
        else
        {
            ExecutorService es;
            List<Callable<Double[][]>> list = new ArrayList<>();
            List<Future<Double[][]>> l2;
            for (int i : LATTICE_SIZES)
            {
                Measure m = new Measure(sweeps, 2.0, i, i);
                list.add(m);
            }

                 es = Executors.newFixedThreadPool(LATTICE_SIZES.length);
                 l2 = es.invokeAll(list);
            try
            {
                es.shutdown();
            }
            catch (Exception e)
            {
                e.printStackTrace();
                System.exit(0);
            }
            for (int i =0;i < LATTICE_SIZES.length; i++)
            {
                try
                {
                    Double[][] o1 = l2.get(i).get();
                    Datapoints.put(LATTICE_SIZES[i], Arrays.deepToString(o1));
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }
        }
        return  Datapoints;
    }

    public static void main(String[] args)
    {
//        double[][] aa = new double[3][3];
//        Arrays.stream(aa).forEach(row -> Arrays.fill(row, new Random().nextBoolean() ? 1:-1));
//        System.out.println(Arrays.deepToString(aa));
        JSONObject o = new JSONObject();
        o.put("temperature", "5.0");
        o.put("sweeps", "1000");
        o.put("row", "100");
        o.put("col", "100");
        o.put("coupling", "1");
        o.put("measure", "100");
        Ising i = new Ising(o);
        double in = 2.0;
        int n = 0;
        Instant s = Instant.now();
        try {
            System.out.println(i.binderCrossing());

        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        System.out.println(Duration.between(s, Instant.now()).toSeconds());
    }


    @Data
    static class Measure implements Callable<Double[][]>
    {
        private int sweeps;
        private double temperature;
        private double binderCumulant;
        private int row;
        private int col;

        private ArrayList<Double[]> binder_values;

        Measure()
        {
            sweeps = 100;
            binderCumulant = 0;
            temperature = 2.0;
        }
        Measure(int s, double t, int r, int c)
        {
            row = r;
            col = c;
            sweeps = s;
            temperature  = t;
            binderCumulant = 0;
            binder_values = new ArrayList<>();

        }
        public Double[][] call()
        {
            Ising i = new Ising();
            i.row = row;
            i.col = col;
            i.lattice = new int[row][col];
            i.sweeps = sweeps;
            i.temperature = temperature;
            i.stack = new Stack<>();
            i.cluster = new boolean[row][col];

            double count = 2.0;
            i.initialiseLattice(row, col);
            while (count < 2.6)
            {
                i.temperature = count;
                i.magnetisation2 = new ArrayList<>();
                i.magnetisation4 = new ArrayList<>();
                i.binderCumulant = 0;
                i.wolff(true);
                binder_values.add(new Double[]{count, i.binderCumulant});


                count+=0.01;
            }
            System.out.println(Arrays.deepToString(binder_values.toArray(new Double[0][0])));
            return binder_values.toArray(new Double[0][0]);

        }
    }
}
