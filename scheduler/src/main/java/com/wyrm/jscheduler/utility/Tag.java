package com.wyrm.jscheduler.utility;

import lombok.Data;

@Data
public class Tag
{

    private String tag;
    private double value;

    public void setTagAndValue(String t, double v)
    {
    tag = t;
    value = v;
    }
    public String toString()
    {
        return "Tag: "+tag+"\tValue: "+value;
    }
}
