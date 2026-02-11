package com.wyrm.jscheduler.utility;


import lombok.Getter;

/**
     * This class is used to represent pairs of data, used in pairwise correlation coefficient calculations
     * @param <T>
     */
    public final class Pair<T, T2>
    {
        @Getter
        private T p1;
        @Getter
        private T2 p2;
        public Pair()
        {}

        public Pair(T p1, T2 p2)
        {
            this.p1 = p1;
            this.p2 = p2;
        }

        public boolean equals(T o1, T2 o2)
        {
            return (p1.equals(o1) && p2.equals(o2)) || (p2.equals(o1) && p1.equals(o2));
        }
        public boolean equals(Pair<T, T2>  o1)
        {
            return (p1.equals(o1.p1) && p2.equals(o1.p2)) || (p2.equals(o1.p1) && p1.equals(o1.p2));
        }



    public String toString()
    {
        return p1+":"+p2;
    }
    }
