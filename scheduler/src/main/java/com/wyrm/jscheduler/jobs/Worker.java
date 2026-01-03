package com.wyrm.jscheduler.jobs;

import com.wyrm.jscheduler.Entity.Output;
import com.wyrm.jscheduler.Entity.Task;
import com.wyrm.jscheduler.repository.OutputRepository;
import com.wyrm.jscheduler.repository.TaskRepository;
import com.wyrm.jscheduler.utility.ANSI;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.json.JSONObject;

import java.sql.SQLException;
import java.time.Instant;
import java.util.Optional;
import java.util.UUID;

@Slf4j
@Data
@NoArgsConstructor
@AllArgsConstructor
public class Worker extends Thread
{
    private Job job;
    private Task task;
    private Output output;
    private TaskRepository taskRepo;
    private OutputRepository outputRepo;

    public void retrieveTask() throws SQLException
    {
        try {
            Optional<Task> optional = taskRepo.findByStatus("PENDING");
            optional.ifPresentOrElse(this::setAndExecute, workerSleep());

        }
        catch (InterruptedException e)
        {
            log.error(ANSI.RED+"Interrupted Exception"+ANSI.RESET);
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

    public void uploadOutput(Output results)
    {

        Output output = new Output(results);
        try
        {
            outputRepo.save(output);
        }
        catch (Exception e)
        {
            log.error(ANSI.RED+"Error for uploading output of task with id: {}, Error: {}"+ANSI.RESET, results.getID(), e.getMessage());
        }
    }
}
