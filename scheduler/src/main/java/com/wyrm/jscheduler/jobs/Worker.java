package com.wyrm.jscheduler.jobs;

import com.wyrm.jscheduler.Entity.Output;
import com.wyrm.jscheduler.Entity.Task;
import com.wyrm.jscheduler.repository.TaskRepository;
import com.wyrm.jscheduler.utility.ANSI;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import java.sql.SQLException;
import java.util.Optional;

@Slf4j
@Data
@NoArgsConstructor
@AllArgsConstructor
public class Worker extends Thread
{
    private Job job;
    private Task task;
    private Output output;
    private TaskRepository service;

    public void retrieveTask() throws SQLException
    {
        try {
            Optional<Task> optional = service.findByStatus("PENDING");
            optional.ifPresentOrElse(this::setAndExecute, workerSleep());

        }
        catch (InterruptedException e)
        {
            log.info(ANSI.RED+"Interrupted Exception"+ANSI.RESET);
        }
    }
    public Runnable workerSleep() throws InterruptedException
    {
        sleep(10000);
        return null;
    }

    public void setAndExecute(Task t)
    {
        task = t;

    }

    public void Execute()
    {

    }
}
