package com.wyrm.jscheduler.repository;
import com.wyrm.jscheduler.Entity.Task;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.rest.core.annotation.RepositoryRestResource;

import java.util.Optional;
import java.util.UUID;

@RepositoryRestResource
public interface TaskRepository extends JpaRepository<Task, UUID>
{
    @Query(value = "SELECT * FROM task WHERE status = ?1 LIMIT 1", nativeQuery = true)
    Optional<Task> findByStatus(String status);
}
