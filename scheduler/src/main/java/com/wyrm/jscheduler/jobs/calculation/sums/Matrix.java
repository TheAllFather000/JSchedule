package com.wyrm.jscheduler.jobs.calculation.sums;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import com.wyrm.jscheduler.jobs.Job;
import com.wyrm.jscheduler.utility.Tag;

@SuppressWarnings("all")
public class Matrix extends RecursiveAction implements Job
{
    public Double[] parallel_matrix;
    private int hi;
    private int low;
    private Double[] matrix_1;


    private Double[] matrix_2;
    int CUTOFF = 100;

    /**
     * Constructor, used primarily for parallel algorithms
     * @param lo
     * @param hi
     * @param o1
     * @param o2
     * @param results
     */
    public Matrix(int lo,  int hi, Double[] o1, Double[] o2, Double[] results)
    {
        this.hi = hi;
        low = lo;
        matrix_1 = (Double[]) o1;
        matrix_2 = (Double[]) o2;
        this.parallel_matrix = results;
    }
    public Matrix(int lo,  int hi, Double[] o1, Double[] results)
    {
        this.hi = hi;
        low = lo;
        matrix_1 = (Double[]) o1;
        this.parallel_matrix = results;
    }

    public Matrix()
    {}

    /**
     * Add together the elements of two single dimension matrices (arrays)
     * @param o Matrix 1
     * @param o_ Matrix 2
     * @return  Singular matrix, representing the sum of the input matrices
     */
    public static Object one_dimensional_sum(Object o, Object o_)
    {
        if (!(o instanceof Number[]) && !(o_ instanceof Number[]))
        {
            return null;
        }
        assert o instanceof Integer[];
        assert o != null;
        Integer[] array = (Integer[]) o;

        assert o_ instanceof Integer[];
        assert o != null;
        Integer[] array2 = (Integer[]) o_;

        if (array.length != array2.length)
        {
            return null;
        }

        Number[] sum = new Number[array.length];
        for (int i = 0; i < array.length; i++)
        {
            int value = array[i] + array2[i];
            sum[i] = value;
        }
        return sum;
    }

    @Override
    public void compute()
    {
        if ((hi - low) < CUTOFF )
        {
            for (int i = low; i < hi; i++)
            {
                parallel_matrix[i] = matrix_1[i] + (matrix_2[i]);
                System.out.println(parallel_matrix[i]);
            }
        }
        else
        {
            Matrix left=  new Matrix(low, (hi+low)/2, matrix_1, matrix_2, parallel_matrix);
            Matrix right = new Matrix((hi+low)/2,hi, matrix_1, matrix_2, parallel_matrix);
            left.fork();
            right.compute();
            left.join();
        }
    }

    /**
     * Add elements of a single matrix together
     * @param o1 2D matrix
     * @return 1D matrix as a result of element-wise addition
     */
    public Double[] elementWise(Double[][] o1)
    {
        assert o1.length !=0;
        if (o1.length == 0)
            return new Double[]{};
        if (o1.length == 1)
        {
            return o1[0];
        }
        for (int  i = 0; i< o1.length;i++)
        {
            if (i == 0)
                continue;
            if (o1[i].length != o1[i-1].length)
                return null;
        }
        Double[] result = new Double[o1[0].length];
        ForkJoinPool fjp = ForkJoinPool.commonPool();
        for (int  i =0; i < o1.length;i++)
        {
            if (i == 0)
                continue;
            Matrix m = new Matrix(0, o1[0].length, o1[i-1], o1[i], result);
            fjp.invoke(m);
        }
        return result;
    }


    /**
     * Add together the contents of two matrices
     * @param o1 Matrix one
     * @param o2 Matrix two
     * @return Input matrices added together
     */
    public Double[][] matrixAddition2D(Double[][] o1, Double[][] o2)
    {
        Double[][] temp = (o1.length != 0) ? o1 : new Double[][]{};

        assert temp!=o1: "Length of matrix 1 is 0";
        temp = (o2.length != 0) ? o2 : new Double[][]{};
        assert temp!=o2: "Length of matrix 2 is 0";
        assert o1.length != o2.length: "Length of matrices are not equal";
        boolean useParallel = false;

        for (int i = 0; i < o1.length; i++)
        {
            assert o1[i].length != o2[i].length;
            if (o1[i].length > 1000) {
                useParallel = true;
            }
        }

        Double[][] results = new Double[o1.length][];
        Double[] temp_array;
        if (useParallel)
        {
            ForkJoinPool fjp = ForkJoinPool.commonPool();
            for (int i = 0; i < o1.length; i++) {
                temp_array = new Double[i];
                Matrix m = new Matrix(0, o1[i].length, o1[i], o2[i], temp_array);
                results[i] = temp_array;
            }
            return results;
        }
        else
        {
            for (int i =0; i < o1.length;i++)
            {
                temp_array = new Double[o1[i].length];
                for (int j = 0; j < o1[i].length;i++)
                {
                    temp_array[i] = o1[i][j] + o2[i][j];
                }
                results[i] = temp_array;
            }
            return results;
        }
    }
    /* valid matrix example
    [1, 2, 3]   * [10, 11]
    [4, 5, 6]     [15, 16]
                  [18, 19]
     */

    /**
     * Performs matrix multiplication on a pair of 2D arrays
     * @param o1 Matrix 1
     * @param o2 Matrix 2
     * @return Matrix Product
     */
    public Double[][] MatrixMultiplication(Double[][] o1, Double[][] o2)
    {
        assert o1.length == 0 || o2.length == 0:"Input matrices are of length 0";
        assert o1.length !=o2[0].length || o1[0].length != o2.length: "Columns on Matrix 1 does not equal rows on matrix 2";
        Double[][] result = new Double[o1.length][o2[0].length];
        for (int row = 0; row < o1.length; row++)
        {
            for (int col = 0; col < o2[0].length; col++)
            {
                result[row][col] = MatrixValue(o1, o2, row, col);
            }
        }
        return result;
    }

    public Double MatrixValue(Double[][] o1, Double[][] o2, int rowVal, int colVal)
    {
        double matrixValue = 0.0;
        Double[] row = o1[rowVal];
        Double[] col = new Double[o2.length];
        for (int m = 0 ; m < o2.length; m++)
        {
            matrixValue+=row[m] + o2[colVal][m];
        }
        return matrixValue;
    }

    public double[][] transpose(double[][] o1)
    {
        assert o1 != null: "Input matrix is null";
        assert o1.length !=0: "Input matrix is empty";
        double[][] o2 = new double[o1[0].length][o1.length];
        for (int rows  = 0 ; rows < o1.length; rows++)
        {
            for (int cols = 0; cols < o1[0].length;cols++)
                o2[cols][rows] = o1[rows][cols];
        }
        return o2;
    }
    public Tag[][] transpose(Tag[][] o1)
    {
        assert o1 != null: "Input matrix is null";
        assert o1.length !=0: "Input matrix is empty";
        Tag[][] o2 = new Tag[o1[0].length][o1.length];
        for (int rows  = 0 ; rows < o1.length; rows++)
        {
            for (int cols = 0; cols < o1[0].length;cols++)
                o2[cols][rows] = o1[rows][cols];
        }
        return o2;
    }


}
