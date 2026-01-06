package com.wyrm.jscheduler.jobs.monte_carlo.financial;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.distribution.NormalDistribution;
import com.wyrm.jscheduler.jobs.Job;
import org.json.JSONObject;

@SuppressWarnings("all")
public class Price implements Job
{
    private Variance v;
    private StandardDeviation sd;
    private Mean m;
    //private ArrayList<Double> returnsList;
    private double[] historical_data;
    private NormalDistribution nd;
    private int repetitions;
    private JSONObject data ;
    public Price(double[] d)
    {
        historical_data = d;
        v = new Variance();
        sd = new StandardDeviation();
        m = new Mean();
        //returnsList = new ArrayList<>();
        nd = new NormalDistribution();
        repetitions = 10;

    }

    public Price(double[] d, int reps)
    {
        historical_data = d;
        v = new Variance();
        sd = new StandardDeviation();
        m = new Mean();
        //returnsList = new ArrayList<>();
        nd = new NormalDistribution();
        repetitions = reps;
    }

    public Price(JSONObject d)
    {
        data = d;
        nd = new NormalDistribution();
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
        assert historical_data.length >1: "Insufficient data for price prediction";
        double[] daily_return = new double[historical_data.length];
        for (int i = 1; i < daily_return.length; i++)
        {
            daily_return[i] = Math.log(historical_data[i]/ historical_data[i-1]);
        }
        v.setBiasCorrected(false);
        sd.setBiasCorrected(false);
        double mean = m.evaluate(daily_return, 0, daily_return.length);
        double variance = v.evaluate(daily_return, mean);
        double drift = mean - (variance/2);
        double rv = Math.sqrt(variance)*nd.inverseCumulativeProbability(Math.random());
        return historical_data[historical_data.length-1]*Math.pow(Math.E, drift+rv);
    }



    /**
     * Returns the predicted price of an item with a user submitted 'drift'.
     * @return predicted price of the asset
     */
    public double pricePrediction(double drift)
    {
        //calculate series of periodic returns
        assert historical_data.length >1: "Insufficient data for price prediction";
        double[] daily_return = new double[historical_data.length];
        for (int i = 1; i < daily_return.length; i++)
        {
            daily_return[i] = Math.log(historical_data[i]/ historical_data[i-1]);
        }
        v.setBiasCorrected(false);
        sd.setBiasCorrected(false);
        double mean = m.evaluate(daily_return, 0, daily_return.length);
        double variance = v.evaluate(daily_return, mean);
        double rv = Math.sqrt(variance)*nd.inverseCumulativeProbability(Math.random());
        return historical_data[historical_data.length-1]*Math.pow(Math.E, drift+rv);
    }


    public double BrowningMotion()
    {

    }
}
