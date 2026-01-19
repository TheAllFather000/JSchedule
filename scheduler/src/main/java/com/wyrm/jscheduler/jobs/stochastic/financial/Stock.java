package com.wyrm.jscheduler.jobs.stochastic.financial;
import lombok.AllArgsConstructor;
import lombok.NoArgsConstructor;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.distribution.NormalDistribution;
import com.wyrm.jscheduler.jobs.Job;
import org.json.JSONObject;
import lombok.Data;
import java.util.Arrays;
import java.util.HashMap;
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
    private double[] returns;
    private double[] innovationSeries;
    private double[] conditionalVariance;

    //used for GARCH prediction
    private double mean;
    private double standardDev;
    private double omega;
    private double alpha;
    private double beta;
    private double logLikelihood;
    private double[] logs;

    //minimum amount of times the MLE is repeated.
    public static final int MIN_DIMENSIONS = 1000;
    public Stock(double[] d)
    {
        historicalData = d;
        v = new Variance();
        sd = new StandardDeviation();
        m = new Mean();
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
        m = new Mean();
        sd = new StandardDeviation();
    }

    public Stock(JSONObject d, double[] returns)
    {
        data = d;
        nd = new NormalDistribution();
        m = new Mean();
        this.returns = returns;
    }

    public Stock(JSONObject d, double[] returns, double mean)
    {
        data = d;
        this.returns = returns;
        this.mean = mean;
    }

    public Stock(JSONObject d, double[] returns, double mean, double o, double a, double b)
    {
        data = d;
        this.returns = returns;
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
            daily_return[i-1] = Math.log(historicalData[i]/ historicalData[i-1]);
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
            daily_return[i-1] = Math.log(historicalData[i]/ historicalData[i-1]);
        }
        v.setBiasCorrected(false);
        sd.setBiasCorrected(false);
        double mean = m.evaluate(daily_return, 0, daily_return.length);
        double variance = v.evaluate(daily_return, mean);
        System.out.println(Arrays.toString(daily_return));
        System.out.println(Arrays.toString(historicalData));
        System.out.println();
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
        returns = new double[historicalData.length-1];
        for (int i = 1; i < historicalData.length; i++) {
            returns[i-1] = Math.log(historicalData[i] / historicalData[i - 1]);
        }
        double drift = m.evaluate(returns);
        v.setBiasCorrected(false);
        double var = v.evaluate(returns);
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
        response.put("logarithm_returns", Arrays.toString(returns));
        response.put("mean_log_returns", drift);
        response.put("predictions", Arrays.toString(predictions));
        return response;
    }


    /**
     *Used for estimating α, ω, β and maximum log-likelihood under Gaussian assumption
     */


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
        returns = new double[historicalData.length-1];
        for (int i = 1; i < historicalData.length; i++) {
            returns[i-1] = Math.log(historicalData[i] / historicalData[i - 1]);
        }
        double drift = m.evaluate(returns);
        v.setBiasCorrected(false);
        double var = v.evaluate(returns);
        for (int i = 0; i < predictions.length; i++)
        {
            predictions[i]  = historicalData[historicalData.length] *
                    Math.pow(Math.E, (drift - var / 2)
                            * ((double) data.get("days")) / 252 + Math.sqrt(var)
                            * Math.sqrt(((double) data.get("days")) / 252));
        }
        response=  new JSONObject();
        response.put("variance", var);
        response.put("logarithm_returns", Arrays.toString(returns));
        response.put("mean_log_returns", drift);
        response.put("predictions", Arrays.toString(predictions));
        return response;
    }


    public JSONObject GARCH()
    {
        assert data != null : "data object is null";
        assert data.has("historical_data") : "Historical data is missing";

        historicalData = Stream.of(((String) data.get("historical_data")).replace("{", "")
                .replace("}", "").replace("[", "")
                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
        returns = new double[historicalData.length - 1];
        innovationSeries = new double[returns.length];
        conditionalVariance = new double[returns.length];
        for (int i = 1; i < historicalData.length; i++)
        {
            returns[i - 1] = Math.log(historicalData[i] / historicalData[i - 1]);
        }
        Processor s = new Processor(data, returns);

        if (data.has("goal"))
        {
            if (((String)data.get("goal")).toLowerCase().equals("minimize"))
            {
                s.setNegativeLogLikelihood(true);
            }
        }
        Thread t = new Thread(s);
        t.start();                                        //MLE is done in another thread.
        double condMean = m.evaluate(returns);         //conditional mean is the average of the log returns
        for (int i = 0; i < innovationSeries.length; i++) {
            innovationSeries[i] = returns[i] - condMean;
        }
        try {
            t.join();
        }
        catch (InterruptedException ie)
        {
            ie.printStackTrace();
        }
        JSONObject output = new JSONObject();
        output.put("alpha", s.getAlpha());
        output.put("beta", s.getBeta());
        output.put("omega", s.getOmega());
        output.put("log_likelihood", s.getLogLikelihood());
        output.put("log_returns", Arrays.toString(returns));
        output.put("innovation_series", Arrays.toString(innovationSeries));
        output.put("conditional_variance", Arrays.toString(s.getCondVariances()));
        output.put("mean_of_log_returns", condMean);
        return output;
    }

    public double[] conditionalVariance(double a, double b, double w)
    {
        conditionalVariance[0] = w/(1-a-b);
        for (int  i = 1; i < conditionalVariance.length;i++)
        {
            conditionalVariance[i] = w+a*Math.pow(returns[i-1],2)+b*conditionalVariance[i-1];
        }
        return conditionalVariance;
    }

    public JSONObject VARHistorical()
    {
        assert data!=null: "JSONObject 'data' is empty";
        historicalData = Stream.of(((String) data.get("historical_data")).replace("{", "")
                .replace("}", "").replace("[", "")
                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
        double[] dailyChanges= new double[historicalData.length-1];
        for (int i = 1 ; i  < dailyChanges.length;i++)
        {
            dailyChanges[i-1] = (historicalData[i]-historicalData[i-1])/historicalData[i-1];
        }
        dailyChanges = quicksort(dailyChanges, 0, dailyChanges.length-1);
        double alpha = (data.has("alpha")) ? (double) data.get("alpha"): 0.05;
        int index = Math.toIntExact(Math.round(alpha * 100));
        System.out.println(index);
        double valueAtRisk = dailyChanges[index];
        double dVar = ((double) data.get("portfolio_value"))*Math.abs(valueAtRisk);
        JSONObject output= new JSONObject();
        output.put("daily_changes", Arrays.toString(dailyChanges));
        output.put("value_at_risk", valueAtRisk);
        output.put("potential_maximum_loss", dVar);
        return output;
    }


    public JSONObject VARVarianceCovariance()
    {
        assert data!=null:"JSONObject 'data' !=null";
        historicalData = Stream.of(((String) data.get("historical_data")).replace("{", "")
                .replace("}", "").replace("[", "")
                .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray();
        double[] dailyChanges= new double[historicalData.length-1];
        for (int i = 1 ; i  < dailyChanges.length;i++)
        {
            dailyChanges[i-1] = (historicalData[i]-historicalData[i-1])/historicalData[i-1];
        }
        Processor s = new Processor(dailyChanges);
        Thread t = new Thread(s);
        t.start();
        double alpha = (data.has("alpha")) ? 1- (double) data.get("alpha"): 0.95;
        try
        {
            t.join();
        }
        catch (InterruptedException ie)
        {
            ie.printStackTrace();
        }
        double valueAtRisk = s.getMean()-s.getStandardDev()*nd.inverseCumulativeProbability(alpha);
        JSONObject output = new JSONObject();
        output.put("returns", Arrays.toString(dailyChanges));
        output.put("mean", s.getMean());
        output.put("standard_deviation", s.getStandardDev());
        output.put("value_at_risk", valueAtRisk);
        return output;
    }

    public JSONObject VARVarianceP()
    {
        assert data!=null:"JSONObject 'data' !=null";
        JSONObject assets = new JSONObject((String) data.get("assets"));
        HashMap<String,double[]> historical = new HashMap<>();
        for (String s: assets.keySet())
        {
            //the asset historical data can exist either as a string or a double[]
            if (assets.get(s).getClass() == double[].class)
                historical.put(s, (double[]) assets.get(s));
            else
                historical.put(s, Stream.of(((String)assets.get(s)).replace("{", "")
                        .replace("}", "").replace("[", "")
                        .replace("]", "").split(",\\s+")).mapToDouble(Double::parseDouble).toArray());
        }
        return null;
    }


    public int partition(double[] a, int l, int h)
    {
        double pivot = a[h];
        int u = l - 1;
        for (int i = l; i <h;i++)
        {
            if (a[i]<=pivot)
            {u++; double temp = a[i]; a[i] = a[u]; a[u] = temp;}

        }
        double temp = a[h]; a[h] = a[u+1];a[u+1] = temp;
        return u+1;
    }
    public double[] quicksort(double[] a, Integer l, Integer h )
    {
        if (h == null)
            h = a.length-1;
        if (l < h)
        {
            int pivot = partition(a, l, h);
            quicksort(a, l,pivot-1 );
            quicksort(a, pivot+1, h);
        }
        return a;
    }

    //values[0]+values[1]*Math.pow((returns[j-1]-mean),2)+values[2]*variances[j-1];
    public static void main(String[] args)
    {
//        Mean m = new Mean();
//        double[][] test = new double[][]{{1,2,3}, {4,5,6}, {7,8,9}};
//        double[] test1 = new double[]{1,2,3,4,5,6,7,8,9};
//        double mean = m.evaluate(test1);
//        double m2 = 0;
//        for (double[] d: test)
//        {
//            m2 +=m.evaluate(d);
//        }
//        m2 /=3;
//        System.out.println(m2 == mean);
        double[] stockPrices = {
                102.45, 98.12, 105.67, 110.34, 107.89, 95.43, 99.87, 103.56,
                108.91, 112.04, 106.78, 101.23, 97.66, 94.88, 96.54, 100.12,
                104.76, 109.33, 113.58, 111.02, 107.41, 103.95, 99.34, 98.76,
                102.88, 105.19, 108.47, 110.92, 114.26, 112.73, 109.05, 104.39,
                101.67, 97.45, 95.92, 98.14, 100.89, 103.21, 106.54, 109.88,
                113.02, 111.47, 108.66, 105.04, 101.98, 99.76, 96.33, 94.57,
                97.89, 100.45, 103.77, 107.12, 110.56, 114.03, 112.18, 109.94,
                106.23, 102.67, 99.98, 97.21, 95.64, 98.32, 101.05, 104.61,
                108.09, 111.83, 115.27, 113.66, 110.48, 107.35, 103.84, 100.29,
                98.55, 96.88, 99.14, 102.73, 105.96, 109.41, 112.95, 116.38,
                114.92, 111.64, 108.72, 105.33, 102.14, 99.67, 97.09, 95.48,
                98.76, 101.38, 104.92, 108.56, 112.07, 115.81, 114.35, 110.96,
                107.84, 104.18, 100.76, 98.21, 96.57, 99.03
        };
        JSONObject o = new JSONObject();
        o.put("historical_data", Arrays.toString(stockPrices));
        o.put("goal","minimize");
        Stock S = new Stock(o);
        o = S.GARCH();
        for (String s: o.keySet())
        {
            System.out.println(s + " : "+o.get(s));
        }
    }

}






@Data
@AllArgsConstructor
class Processor implements Runnable
{
    private double mean;
    private double standardDev;
    private double omega;
    private double alpha;
    private double beta;
    private double logLikelihood;
    private double[] logs;
    private double[] returns;
    private double[] asset_mean;
    private boolean negativeLogLikelihood;
    private double[] condVariances;


    private JSONObject data;
    private HashMap<String,double[]> assetReturns;
    private Mean m;
    private StandardDeviation sd;
    private NormalDistribution nd;

    public Processor()
    {
        m = new Mean();
        sd = new StandardDeviation();
        nd = new NormalDistribution();

    }

    public Processor(double[] r)
    {
        m = new Mean();
        sd = new StandardDeviation();
        nd = new NormalDistribution();
        returns = r;
    }
    public Processor(JSONObject d)
    {
        m = new Mean();
        sd = new StandardDeviation();
        nd = new NormalDistribution();
        data = d;
    }
    public Processor(JSONObject d, double[] r)
    {
        m = new Mean();
        sd = new StandardDeviation();
        nd = new NormalDistribution();
        data = d;
        returns = r;
    }
    @Override
    public void run()
    {
        if (returns !=null && data == null)
        {
            mean = m.evaluate(returns);
            standardDev = sd.evaluate(returns);
            return;
        }

        else if (assetReturns !=null)
        {
            asset_mean = new double[assetReturns.size()];
            int i = 0;
            for (String S: assetReturns.keySet())
            {
                asset_mean[i] = m.evaluate(assetReturns.get(S));
                i++;
            }
            return;
        }

        /**
         *Used for estimating α, ω, β and maximum log-likelihood under Gaussian assumption
         */
        else if (data !=null && returns !=null)

        {
            double w = 0.00001;
            double a = 0.05;
            double b = 0.85;
            double totalLogs = 0;
            MLE mle = new MLE(returns, mean);
            SimplexOptimizer op = new SimplexOptimizer(new SimpleValueChecker(1e-12, 1e-12));
            NelderMeadSimplex nd = new NelderMeadSimplex(new double[]{1e-9, 0.02, 0.05});
            PointValuePair pv;
            if (negativeLogLikelihood)
            {
                mle.setLogSign(true);
                if (data.has("dimensions")) pv = op.optimize(new MaxEval((int) data.get("dimensions")), new ObjectiveFunction(mle), new InitialGuess(new double[]{1e-6, 0.10, 0.90}), GoalType.MINIMIZE, nd);
                else pv = op.optimize(new MaxEval(Stock.MIN_DIMENSIONS), new ObjectiveFunction(mle), new InitialGuess(new double[]{1e-6, 0.10, 0.90}), GoalType.MINIMIZE, nd);
                double[] p = pv.getPoint();
                omega = p[0];
                alpha = p[1];
                beta = p[2];
                logLikelihood = -pv.getValue();
                condVariances = mle.getConditionalVariance();
            }
            else
            {
                if (data.has("dimensions")) pv = op.optimize(new MaxEval((int) data.get("dimensions")), new ObjectiveFunction(mle), new InitialGuess(new double[]{1e-6, 0.10, 0.90}), GoalType.MAXIMIZE, nd);
                else pv = op.optimize(new MaxEval(Stock.MIN_DIMENSIONS), new ObjectiveFunction(mle), new InitialGuess(new double[]{1e-6, 0.10, 0.90}), GoalType.MAXIMIZE, nd);
                double[] p = pv.getPoint();
                omega = p[0];
                alpha = p[1];
                beta = p[2];
                logLikelihood = pv.getValue();
                condVariances = mle.getConditionalVariance();
            }
        }

    }
    public void setSign(boolean b)
    {
        negativeLogLikelihood = b;
    }

}
