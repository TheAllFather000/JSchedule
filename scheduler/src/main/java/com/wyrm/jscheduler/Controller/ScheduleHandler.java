package com.wyrm.jscheduler.Controller;
import com.wyrm.jscheduler.Entity.Task;
import org.json.JSONObject;

import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import com.wyrm.jscheduler.Service.TaskService;
import java.time.Instant;
import java.util.UUID;

@Slf4j
@RestController
@RequestMapping("/api/task")
@CrossOrigin()
public class ScheduleHandler
{
    @Autowired
    private TaskService taskService;

    @PostMapping("/submit")
    public JSONObject submitTask(@RequestBody Task task)
    {
        UUID id = UUID.randomUUID();
        task.setID(id);
        task.setTimestamp(Instant.now());
        taskService.save(task);
        JSONObject response = new JSONObject();
        response.put("task_id", id);
        response.put("status", "PENDING");
        response.put("job_type", task.getType());
        response.put("data", task.getData());
        response.put("timestamp", task.getTimestamp());
        return response;
    }

    @GetMapping("/status")
    public JSONObject checkTaskStatus(@RequestBody UUID id)
    {
        Task t = taskService.getReferenceById(id);
        JSONObject response = new JSONObject();
        response.put("task_id", t.getID());
        response.put("job_type",t.getType());
        response.put("status", t.getStatus());
        response.put("submitted_at", t.getTimestamp());
        return response;
    }
}
