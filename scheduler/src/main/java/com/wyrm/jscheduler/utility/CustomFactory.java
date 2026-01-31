import java.util.concurrent.ThreadFactory;
public class CustomFactory implements ThreadFactory
{
    private final AtomicInteger count = new AtomicInteger(0);

    @Override
    public Thread newThread(Runnable r)
    {
        Thread t = new Thread(r);
        t.setName("Ising: " +count.get())
        count.addAndGet(1);
        return t;
    }
}
