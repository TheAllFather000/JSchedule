package com.wyrm.jscheduler.utility;

import java.util.*;




public class Pairs<T>
{
    private final ArrayList<Pair<T>> pairs;
    public Pairs(){pairs = new ArrayList<>(); }


    public List<Pair<T>> findPairs(Set<T> keys)
    {
        Arrays.sort(keys.toArray());
            boolean present = false;
            List<Pair<T>> pairs_ = new ArrayList<>();
            for (T t: keys)
            {
                for (T tt: keys)
                {
                    if (!t.equals(tt))
                    {
                        Pair<T> p1;
                        if (((String) t).compareTo(((String) tt)) <0)
                            p1 = new Pair<>(t, tt);
                        else
                            p1 = new Pair<>(tt, t);
                        for (Pair<T> p : pairs_)
                        {
                            if (p.equals(t,tt))
                            {
                                present = true;
                            }
                        }
                        if (!present) pairs_.add(p1);
                    }
                    present = false;
                }
            }
        System.out.println(pairs_);
        return pairs_;
    }

    public static void main(String[] args)
    {
        for (int  i = 0 ; i < 10; i++)
            for (int j = 0 ; j < 10;j++)
                if (i == j)
                    System.out.println(i + " "+ j);
    }

}