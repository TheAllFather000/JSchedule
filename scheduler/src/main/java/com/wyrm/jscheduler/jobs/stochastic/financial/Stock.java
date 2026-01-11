package com.wyrm.jscheduler.jobs.stochastic.financial;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.distribution.NormalDistribution;
import com.wyrm.jscheduler.jobs.Job;
import org.json.JSONObject;
import net.finmath.timeseries.models.parametric.GARCH;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

@SuppressWarnings("all")
public class Stock extends Thread implements Job
{
    private Mean m;
    private Variance v;
    private StandardDeviation sd;
    private NormalDistribution nd;

    private JSONObject data ;
    private int repetitions;
    private double[] historicalData;
    private double[] logReturns;
    private double[] innovationSeries;

    //used for GARCH prediction
    private double mean;
    private double omega;
    private double alpha;
    private double beta;
    private double logTolerance;
    private double[] logs;
    private GARCH garch;

    //minimum amount of times the MLE is repeated.
    public static final int MIN_DIMENSIONS = 100;
    public Stock(double[] d)
    {
        historicalData = d;
        v = new Variance();
        sd = new StandardDeviation();
        m = new Mean();
        //returnsList = new ArrayList<>();
        nd = new NormalDistribution();
        repetitions = 10;

    }

    public Stock(double[] d, int reps)
    {
        historicalData = d;
        v = new Variance();
        sd = new StandardDeviation();
        m = new Mean();
        //returnsList = new ArrayList<>();
        nd = new NormalDistribution();
        repetitions = reps;
    }

    public Stock(JSONObject d)
    {
        data = d;
        nd = new NormalDistribution();
    }

    public Stock(JSONObject d, int reps)
    {
        data = d;
        repetitions =  reps;
        nd = new NormalDistribution();
        m = new Mean();
    }
    public Stock(GARCH g)
    {
        garch = g;
    }

    public Stock(JSONObject d, double[] logReturns, double mean)
    {
        data = d;
        this.logReturns = logReturns;
        nd = new NormalDistribution();
        m = new Mean();
        this.mean = mean;
    }

    public Stock(JSONObject d, double[] logReturns, double mean, double o, double a, double b)
    {
        data = d;
        this.logReturns = logReturns;
        nd = new NormalDistribution();
        m = new Mean();
        this.mean = mean;
        omega = o;
        alpha = a;
        beta = b;
    }

    public double[] RepeatPricePrediction()
    {
        double[] results = new double[repetitions];
        for (int i = 0; i < repetitions; i++)
        {
            results[i] = pricePrediction();
        }
        return results;
    }


    /**
     * Returns a range of price predictions using a specific 'drift' value
     * @param drift
     * @return price predictions
     */
    public double[] RepeatPricePrediction(double drift)
    {
        double[] results = new double[repetitions];
        for (int i = 0; i < repetitions; i++)
        {
            results[i] = pricePrediction(drift);
        }
        return results;
    }


    /**
     * Returns the predicted price of an item using historical data of an asset
     * @return predicted price of the asset
     */
    public double pricePrediction()
    {
        //calculate series of periodic returns
        assert historicalData.length >1: "Insufficient data for price prediction";
        double[] daily_return = new double[historicalData.length-1];
        for (int i = 1; i < daily_return.length; i++)
        {
            daily_return[i] = Math.log(historicalData[i]/ historicalData[i-1]);
        }
        v.setBiasCorrected(false);
        sd.setBiasCorrected(false);
        double mean = m.evaluate(daily_return, 0, daily_return.length);
        double variance = v.evaluate(daily_return, mean);
        double drift = mean - (variance/2);
        double rv = Math.sqrt(variance)*nd.inverseCumulativeProbability(Math.random());
        return historicalData[historicalData.length-1]*Math.pow(Math.E, drift+rv);
    }

    /**
     * Returns the predicted price of an item with a user submitted 'drift'.
     * @return predicted price of the asset
     */
    public double pricePrediction(double drift)
    {
        //calculate series of periodic returns
        assert historicalData.length >1: "Insufficient data for price prediction";
        double[] daily_return = new double[historicalData.length-1];
        for (int i = 1; i < daily_return.length; i++)
        {
            daily_return[i] = Math.log(historicalData[i]/ historicalData[i-1]);
        }
        v.setBiasCorrected(false);
        sd.setBiasCorrected(false);
        double mean = m.evaluate(daily_return, 0, daily_return.length);
        double variance = v.evaluate(daily_return, mean);
        double rv = Math.sqrt(variance)*nd.inverseCumulativeProbability(Math.random());
        return historicalData[historicalData.length-1]*Math.pow(Math.E, drift+rv);
    }


    public JSONObject GBM()
    {
        double[] predictions;
        JSONObject response;
        assert data != null: "JSONObject 'data' is empty";
        assert data.has("stock_prices"): "Historical data is missing.";
        if (data.has("volumes"))
            predictions = new double[(int) data.get("volumes")];
        else
            predictions = new double[10];
        historicalData = Stream.of(((String) data.get("historical")).replace("{", "")
                .replace("}", "").replace("[", "")
                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
        logReturns = new double[historicalData.length-1];
        for (int i = 1; i < historicalData.length; i++) {
            logReturns[i] = Math.log(historicalData[i] / historicalData[i - 1]);
        }
        double drift = m.evaluate(logReturns);
        v.setBiasCorrected(false);
        double var = v.evaluate(logReturns);
        for (int i = 0; i < predictions.length; i++)
        {
            predictions[i]  = historicalData[historicalData.length] *
                    Math.pow(Math.E, (drift - var / 2)
                            * ((double) data.get("days")) / 252 + Math.sqrt(var)
                            * Math.sqrt(((double) data.get("days")) / 252)
                            * nd.inverseCumulativeProbability(Math.random()));

        }
        response=  new JSONObject();
        response.put("variance", var);
        response.put("logarithm_returns", Arrays.toString(logReturns));
        response.put("mean_log_returns", drift);
        response.put("predictions", Arrays.toString(predictions));
        return response;
    }

    //used for estimating α, ω, β and log-likelihood under Gaussian assumption
    public void start()
    {
        //initial values
        double w =  0.00001;
        double a = 0.05;
        double b = 0.85;
        double totalLogs = 0;
        double[] condVar = new double[logReturns.length];
        /*the values of omega, alpha and beta are randomly generated within a restricted range:
         α∈[0.05,0.15], β∈[0.85,0.95]
         */
        Random r = new Random();
        logs = new double[logReturns.length];

        condVar[0] = w/(1-a-b);
        for (int i = 1; i < condVar.length; i++)
        {
            condVar[i] = w+a*Math.pow(logReturns[i],2) + b*condVar[i-1];
            logs[i] = (-1/2)*(Math.log(2*Math.PI)+Math.log(condVar[i])+(logReturns[i]/condVar[i]));
            totalLogs +=logs[i];
        }
        double previousTotal = totalLogs;
        if (data.has("dimensions"))
        {
            int d = 0;
            while (d < (int) data.get("dimensions"))
            {
                a = r.nextDouble(0.05,0.20);
                b = r.nextDouble(0.80,0.95);
                if (a+b >=1)
                    a = alpha_beta(a,b,r);

            }
        }

    }

    public double alpha_beta(double a, double b, Random r)
    {

            while ((a + b) >=1)
                a = r.nextDouble(0.05,0.20);
            return a;
    }

    public JSONObject GBMRolling()
    {

        assert data != null: "data object is null";
        assert data.has("historical_data"): "Historical data is missing.";
        double[] predictions;
        JSONObject response;
        if (data.has("data"))
            predictions = new double[(int) data.get("days")];
        else
            predictions = new double[10];
        historicalData = Stream.of(((String) data.get("historical")).replace("{", "")
                .replace("}", "").replace("[", "")
                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
        logReturns = new double[historicalData.length-1];
        for (int i = 1; i < historicalData.length; i++) {
            logReturns[i-1] = Math.log(historicalData[i] / historicalData[i - 1]);
        }
        double drift = m.evaluate(logReturns);
        v.setBiasCorrected(false);
        double var = v.evaluate(logReturns);
        for (int i = 0; i < predictions.length; i++)
        {
            predictions[i]  = historicalData[historicalData.length] *
                    Math.pow(Math.E, (drift - var / 2)
                            * ((double) data.get("days")) / 252 + Math.sqrt(var)
                            * Math.sqrt(((double) data.get("days")) / 252)
                            * nd.inverseCumulativeProbability(Math.random()));

        }
        response=  new JSONObject();
        response.put("variance", var);
        response.put("logarithm_returns", Arrays.toString(logReturns));
        response.put("mean_log_returns", drift);
        response.put("predictions", Arrays.toString(predictions));
        return response;
    }


//    public JSONObject GARCH()
//    {
//        assert data != null: "data object is null";
//        assert data.has("historical_data"): "Historical data is missing";
//
//        historicalData = Stream.of(((String) data.get("historical")).replace("{", "")
//                .replace("}", "").replace("[", "")
//                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
//        logReturns = new double[historicalData.length-1];
//        innovationSeries = new double[logReturns.length];
//        for (int i = 1; i < historicalData.length; i++)
//        {
//            logReturns[i-1] = Math.log(historicalData[i] / historicalData[i - 1]);
//        }
//        //conditional mean is the average of the log returns
//        double condMean = m.evaluate(logReturns);
//        for (int i = 0; i < innovationSeries.length; i++)
//        {
//            innovationSeries[i] = logReturns[i] - condMean;
//        }
//
//    }


    public static void main(String[] args)
    {
        double[] s = new double[100];
        for (int  i = 0; i < s.length;i++)
        {
            if ( i %2 ==0)
               s[i]= 50 + (150 + 50) * new Random().nextDouble();
            else
                s[i] = 50 + (150 - 50) * new Random().nextDouble();
        }
        GARCH g = new GARCH(s);
        Map<String, Object> g2 = new HashMap<>();
        g2.put("Omega", 0.00001);
        g2.put("Alpha", 0.15);
        g2.put("Beta", 0.80);
        Map<String, Object> g1 = g.getBestParameters();
        Map<String, Object> g4 = g.getBestParameters(g2);
        Map<String, Object> g3 = g.getBestParameters();
        System.out.println(g1.toString());
        System.out.println(g4.toString());
        System.out.println(g3.toString());
        Random r= new Random();

        for (int i = 0; i < 100; i++)
        {
            System.out.println(r.nextDouble(0.05, 0.20));
        }
    }
}
