package com.wyrm.jscheduler.repository;

import com.wyrm.jscheduler.Entity.Output;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.rest.core.annotation.RepositoryRestResource;

import java.util.UUID;

@RepositoryRestResource
public interface OutputRepository extends JpaRepository<Output, UUID>
{

}
