package com.wyrm.jscheduler.Controller;
import com.wyrm.jscheduler.Entity.Output;
import com.wyrm.jscheduler.Entity.Task;
import com.wyrm.jscheduler.repository.OutputRepository;
import org.json.JSONObject;
import com.wyrm.jscheduler.utility.ANSI;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import com.wyrm.jscheduler.repository.TaskRepository;
import java.time.Instant;
import java.util.Optional;
import java.util.UUID;

@Slf4j
@RestController
@RequestMapping("/api/task")
@CrossOrigin()
public class ScheduleHandler
{
    @Autowired
    private TaskRepository taskService;

    @Autowired
    private OutputRepository outputService;

    @PostMapping("/submit")
    public JSONObject submitTask(JSONObject object) {
        JSONObject response = new JSONObject();
        UUID id = UUID.randomUUID();
        try {
            Task task = new Task();
            task.setID(id);
            if (!object.keySet().contains("data")) {
                response.put("Key Error", "key-value pair for 'data' is missing");
                return response;
            }
            if (!object.keySet().contains("type")) {
                response.put("Key Error", "key-value pair for 'task' is missing");
                return response;
            }
            task.setData((JSONObject) object.get("data"));
            task.setTask((String) object.get("type"));
            task.setTimestamp(Instant.now());
            taskService.save(task);
            response.put("id", id);
            response.put("status", "PENDING");
            response.put("task", task.getTask());
            response.put("data", task.getData());
            response.put("timestamp", task.getTimestamp());
            log.info(ANSI.BLUE+"Submitted task with ID: {}"+ANSI.RESET, id);
            return response;
        }
        catch (Exception e)
        {
            log.error(ANSI.RED+"Error for submitting task with id: {}, Error: {}"+ANSI.RESET,id, e.getMessage());
            return response;
        }
    }

    @GetMapping("/status")
    public String checkTaskStatus(@RequestParam("task_id") String id)
    {
        JSONObject response = new JSONObject();
        try {
            Optional<Task> op = taskService.findById(UUID.fromString(id));
            if (op.isPresent()) {
                Task t = op.get();
                response.put("task_id", t.getID());
                response.put("job_type", t.getTask());
                response.put("status", t.getStatus());
                response.put("submitted_at", t.getTimestamp());
                log.info(ANSI.BLUE+"Status request for task with id {}"+ANSI.RESET, id);
            }
            else
            {
                response.put("Error", "No task with id:{"+id+'}');
                log.info(ANSI.MAGENTA+"No task with id {}"+ANSI.RESET, id);

            }
            return response.toString();
        }
        catch (Exception e)
        {
            log.error(ANSI.RED+"Error for status request of task with id: {}, Error: {}"+ANSI.RESET, id, e.getMessage());
            return response.toString();
        }
    }

    @GetMapping("/retrieve")
    public String retrieveOutput(@RequestParam("task_id") String id)
    {
        JSONObject response = new JSONObject();
        try
        {
            Optional<Output> op = outputService.findById(UUID.fromString(id));
            if (op.isPresent())
            {
                Output out = op.get();
                response.put("task_id", out.getID());
                response.put("attempts", out.getAttempts());
                response.put("results", out.getResults());
                response.put("processing_time", out.getProcessing_time());
                response.put("execution_time", out.getExecution_time());
                response.put("execution_time", out.getCompletion_time());
            }
            else
            {
                response.put("Error", "No task with id:{"+id+'}');
                log.info(ANSI.BLUE+"Output request for task with id {}"+ANSI.RESET, id);

            }
            return response.toString();
        }
        catch (Exception e)
        {
            log.error(ANSI.RED+"Error for output request of task with id: {}, Error: {}"+ANSI.RESET, id, e.getMessage());
            return response.toString();
        }
    }
}
