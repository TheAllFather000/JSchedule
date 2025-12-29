package com.wyrm.jscheduler.Service;
import com.wyrm.jscheduler.Entity.Task;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.rest.core.annotation.RepositoryRestResource;
import org.springframework.stereotype.Service;

import java.util.UUID;

@RepositoryRestResource
public interface TaskService extends JpaRepository<Task, UUID>
{

}
