package com.wyrm.jscheduler.utility;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
public class CustomFactory implements ThreadFactory
{
    private final AtomicInteger count = new AtomicInteger(0);

    private boolean daemon;
    public CustomFactory()
    {

    }
    public CustomFactory(boolean daemon)
    {
        this.daemon = daemon;
    }
    @Override
    public Thread newThread(Runnable r)
    {
        Thread t = new Thread(r);
        t.setDaemon(daemon);
        t.setName("Ising: " +count.get());
        count.addAndGet(1);
        return t;
    }
}
