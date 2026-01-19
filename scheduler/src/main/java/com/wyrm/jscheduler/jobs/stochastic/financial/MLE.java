package com.wyrm.jscheduler.jobs.stochastic.financial;

import org.apache.commons.math3.analysis.MultivariateFunction;

@SuppressWarnings("all")
public class MLE implements MultivariateFunction
{
    private double[] returns;
    private double mean;
    private double[] variances;
    private double[] condVariances;
    private boolean returnNegativeLog = false;
    public MLE (double[] returns, double mean)
    {
        this.returns  = returns;
        this.mean = mean;
    }

    //used for minimising of log-likelihood
    @Override
    public double value(double[] values)
    {
        if (values[0] <=0 || values[1] <0  || values[2] < 0 || values[1]+values[2] >=1) return 1e100;
        variances = new double[returns.length];
        variances[0] = values[0]/(1-values[1]-values[2]);
        if (variances[0] <=0 || Double.isNaN(variances[0])) return 1e100;
        double logs = 0;
        for (int j = 1; j < variances.length;j++)
        {
            variances[j] = values[0]+values[1]*Math.pow((returns[j-1]-mean),2)+values[2]*variances[j-1];
            if (variances[j] <=0 || Double.isNaN(variances[j])) return 1e100;
            double current  = returns[j] - mean;
            logs+= -0.5*(Math.log(2.0*Math.PI)+Math.log(variances[j])+(current*current)/variances[j]);
        }
        return returnNegativeLog ? -logs:logs;
    }


    public void setLogSign(boolean t)
    {
        returnNegativeLog = t;
    }

    public double[] getConditionalVariance()
    {return variances;}
}
