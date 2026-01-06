package com.wyrm.jscheduler.jobs.monte_carlo.math;

import java.util.Random;
import com.wyrm.jscheduler.jobs.Job;
import lombok.Data;

@Data
public class PI implements Job
{
    private int circle_points = 0;
    private int square_points = 0;
    private int interval = 0;
    private final Random r = new Random();
    public static final int ITERATIONS = 10000;

    public double EstimatePi() {
        while (true)
        {
            if (interval >= ITERATIONS*ITERATIONS)
                return 4 *  ( (double) circle_points / square_points);
            double x_point = Math.random()*2-1;
            double y_point = Math.random()*2-1;
            double d = x_point * x_point + y_point * y_point;
            if (d <= 1) circle_points++;
            square_points++;
            interval++;
        }

    }



}
